#!/usr/bin/env python3
"""Mission Creek test cases for the OBWB Raven reach-loss patch (decision D181; method doc
docs/OBWB_LosingReachMethod_v1_20261031.md Sections 6 and 9 in mission-creek-mass-balance).

Usage:
  python3 reachloss_mc.py build      # write the case model folders (never touches existing models)
  python3 reachloss_mc.py run        # run every case (reference binary for D0ref, TEST binary for the rest)
  python3 reachloss_mc.py qaqc       # first-principles and reproduction checks, one PASS/FLAG line each; exits 1 on FLAG

Acceptance (stated before any run, method doc 6.4 and 9.2):
  A1  D0 (v2, TEST binary) byte-identical to D0ref (v2, unpatched reference binary) in every output file.
  A2  D1 (D32 lookup table replaced by :ReachLoss EMPIRICAL_Q 0.17 1 0.55, SINK) reproduces v2 Hydrographs
      within print precision (max |diff| <= 1e-4 m3/s, all columns) and the D32 law day by day (1e-9).
      [E58: the first version used F80's 3.271 GL/yr, which is from bdr_wm_dem, not v2.]
  A3  D2 (HRUGROUP SOIL[2] of the SB209 HRUs) and D3 (SUBBASIN 201, k 56 d): end-of-run watershed MB error differs
      from D1's by <= 1 mm (isolates the new booking; v2's raw MB error carries the D113 spill term).
  A4  Supply limit (P4): Qact <= pre-diversion flow every day; no negative flow anywhere (P4).
  A5  Volume two ways (P1, P12): ReachLoss.csv vol_step = 0.5 (Qact + Qact_prev) 86400 every day (1e-6 relative).
  A6  D3 return store obeys the linear-store law (P6): store(t+1) = store(t) exp(-1/k) + vol_step(t+1) within 1e-9 relative.
  A7  Loss law (P2, P3) on the short cases D4 (DISCONNECTED), D5a (CONNECTED, water table below the bed bottom),
      D5b (CONNECTED, h_gw 1 m above the bed): Qact(t) equals the Python reference law evaluated at the start-of-step
      flow Qpre written by the patch within 1e-8 m3/s; D5a equals D4 exactly (disconnection limit, P3); D5b loss is 0 whenever
      stage is below h_gw and gain_clip_steps > 0 (no flow against the gradient, P2).
  P8 (no step at the season boundary) is reported for D1 as FLAG KNOWN: D1 must reproduce D32, which has that step.
Test parameters for D4/D5 (Kbed 0.214 m/d from F202 kbed_MP; d 1 m; length 3000 m; W = 8 Q^0.5, h = 0.3 Q^0.4;
z_bed 400 m) are test values for the code checks only. They are not calibrated and produce no results for use.
"""
import os, sys, shutil, subprocess, filecmp, math, glob
import pandas as pd
import numpy as np

MC = os.path.expanduser("~/Projects/mission-creek-mass-balance")
BASE = os.path.join(MC, "models/v2_working")   # D151 snapshot of bdr_wm2_dem_v412 (real files)
FORCING = os.path.join(MC, "models/cal2_x108_2024_tqm_v412")   # forcing NetCDFs named in the v2 manifest
STEM = "natural-calibrations-mission-nat-v1-1-0-cal1-validation"
REF_BIN = os.path.expanduser("~/Projects/raven-macos/dist/raven-v4.12-macos-arm64/Raven")
TEST_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dist", "raven-v4.12-reachloss-TEST-macos-arm64", "Raven")
MODELS = os.path.join(MC, "models")
RUNS = os.path.join(MC, "outputs/runs")
P = "rl181"
CASES = ["D0ref", "D0", "D1", "D2", "D3", "D4", "D5a", "D5b"]
SHORT = {"D4", "D5a", "D5b"}
SHORT_END = "1996-12-31 00:00:00"
RL_PAR_DISC = [0.214, 1.0, 3000.0, 8.0, 0.5, 0.3, 0.4, 0.0]

def rl_block(case):
    if case == "D1":
        return ":ReachLoss 209 EMPIRICAL_Q 166 242\n  :Parameters 0.17 1.0 0.55\n  :Destination SINK\n:EndReachLoss\n"
    if case == "D2":
        return ":ReachLoss 209 EMPIRICAL_Q 166 242\n  :Parameters 0.17 1.0 0.55\n  :Destination HRUGROUP RL_Fan209 SOIL[2]\n:EndReachLoss\n"
    if case == "D3":
        return ":ReachLoss 209 EMPIRICAL_Q 166 242\n  :Parameters 0.17 1.0 0.55\n  :Destination SUBBASIN 201 56\n:EndReachLoss\n"
    par = " ".join(str(x) for x in RL_PAR_DISC)
    if case == "D4":
        return f":ReachLoss 209 DISCONNECTED 166 242\n  :Parameters {par}\n  :Destination SINK\n:EndReachLoss\n"
    if case == "D5a":
        return f":ReachLoss 209 CONNECTED 166 242\n  :Parameters {par} 400.0 390.0\n  :Destination SINK\n:EndReachLoss\n"
    if case == "D5b":
        return f":ReachLoss 209 CONNECTED 166 242\n  :Parameters {par} 400.0 401.0\n  :Destination SINK\n:EndReachLoss\n"
    return None

def law(case, Q):
    """Python reference of CalcReachLossRate (independent re-implementation, same equations as the spec)."""
    Q = np.maximum(Q, 0.0)
    K, b, ln, aW, bW, aH, bH, cap = RL_PAR_DISC
    W = aW * Q ** bW; h = aH * Q ** bH
    if case == "D4":
        dh = h + b
    else:
        zb, hg = 400.0, (390.0 if case == "D5a" else 401.0)
        dh = (zb + h) - np.maximum(hg, zb - b)
    L = K * dh / b * W * ln / 86400.0
    L = np.where(L < 0, 0.0, L)
    L = np.where(Q <= 0, 0.0, L)
    return np.minimum(L, Q)

def qa(label, ok, flags, known=False):
    tag = "PASS" if ok else ("FLAG KNOWN" if known else "FLAG")
    print(f"{label} [{tag}]")
    if not ok and not known:
        flags.append(label)

def build():
    flags = []
    real = {f: os.path.realpath(os.path.join(BASE, f)) for f in os.listdir(BASE)}
    for nc in ("pr", "tasmax", "tasmin"):
        real[f"{nc}.HRU.timeseries.V1.0.1.nc"] = os.path.realpath(os.path.join(FORCING, f"{nc}.HRU.timeseries.V1.0.1.nc"))
    missing = [f for f, r in real.items() if not os.path.exists(r)]
    qa(f"every model file resolves to an existing file (missing: {missing})", not missing, flags)
    rvh = open(real[f"{STEM}.rvh"]).read()
    lines = rvh.split("\n")
    i0 = next(i for i, l in enumerate(lines) if l.startswith(":HRUs"))
    i1 = next(i for i, l in enumerate(lines) if l.startswith(":EndHRUs"))
    ids = [l.split(",")[0].strip() for l in lines[i0 + 3:i1] if l.strip() and not l.strip().startswith("#") and l.split(",")[5].strip() == "209"]
    qa(f"SB209 HRUs found for the RL_Fan209 group: {len(ids)}", len(ids) > 0, flags)
    mo = open(real["managed_operations.rvt"]).read()
    a = mo.find(":FlowDiversionLookupTable 209 -1 166 242"); e = mo.find(":EndFlowDiversionLookupTable", a)
    qa("exactly one D32 lookup table in managed_operations.rvt", mo.count(":FlowDiversionLookupTable 209 -1") == 1 and a >= 0 and e > a, flags)
    if flags: sys.exit("build stopped: " + "; ".join(flags))
    e += len(":EndFlowDiversionLookupTable")
    for c in CASES:
        d = os.path.join(MODELS, f"{P}_{c}_v412")
        if os.path.exists(d): sys.exit(f"{d} exists; not overwriting")
        os.makedirs(d)
        for f, tgt in real.items():
            os.symlink(tgt, os.path.join(d, f))
        blk = rl_block(c)
        if blk:
            os.remove(os.path.join(d, "managed_operations.rvt"))
            open(os.path.join(d, "managed_operations.rvt"), "w").write(mo[:a] + "# OBWB D181 test: D32 lookup table replaced by :ReachLoss\n" + blk + mo[e:])
        if c == "D2" or c in SHORT:
            os.remove(os.path.join(d, f"{STEM}.rvi"))
            rvi = open(real[f"{STEM}.rvi"]).read()
            if c == "D2":
                rvi += "\n:DefineHRUGroups RL_Fan209\n"
            if c in SHORT:
                rvi = "\n".join((":EndDate\t" + SHORT_END) if l.startswith(":EndDate") else l for l in rvi.split("\n"))
            open(os.path.join(d, f"{STEM}.rvi"), "w").write(rvi)
        if c == "D2":
            os.remove(os.path.join(d, f"{STEM}.rvh"))
            grp = ":HRUGroup RL_Fan209\n" + ",".join(ids) + "\n:EndHRUGroup\n"
            open(os.path.join(d, f"{STEM}.rvh"), "w").write(rvh + "\n" + grp)
        print(f"built {d}")
    print("build: 0 FLAG")

TEST_SHA = "6ff517369ca942c2759068b7f2f29029bf9f5ec8ca351e7e92b7c73a5e50e192"   # D181 hardened build; see BUILD_INFO.txt
REF_SHA = "b95accb90e3f76fafb77ce705b5b3077e1025f09bc09147af92a8c27ea961b4d"    # unpatched reference build (raven-macos)

def sha256(path):
    import hashlib
    return hashlib.sha256(open(path, "rb").read()).hexdigest()

def run():
    flags = []
    qa("TEST binary sha256 is the recorded hardened build", sha256(TEST_BIN) == TEST_SHA, flags)
    qa("reference binary sha256 is the recorded unpatched build", sha256(REF_BIN) == REF_SHA, flags)
    if flags: sys.exit("run stopped: " + "; ".join(flags))
    procs = []
    for c in CASES:
        out = os.path.join(RUNS, f"{P}_{c}")
        if os.path.exists(out): print(f"{out} exists; skipping"); continue
        os.makedirs(out)
        b = REF_BIN if c == "D0ref" else os.path.abspath(TEST_BIN)
        log = open(os.path.join(out, "raven_stdout.txt"), "w")
        procs.append((c, subprocess.Popen(["caffeinate", "-i", "-s", b, STEM, "-o", out + "/"], cwd=os.path.join(MODELS, f"{P}_{c}_v412"), stdout=log, stderr=subprocess.STDOUT)))
        while sum(p.poll() is None for _, p in procs) >= 4:
            import time; time.sleep(10)
    for c, p in procs:
        p.wait(); print(f"{c} exit {p.returncode}")

def rd(c, name):
    d = pd.read_csv(os.path.join(RUNS, f"{P}_{c}", f"{STEM}_{name}.csv"), dtype=str)
    for col in d.columns:
        if col not in ("date", "hour"):
            d[col] = pd.to_numeric(d[col].replace("---", np.nan))
    return d

def mb_end(c):
    w = rd(c, "WatershedStorage")
    col = [x for x in w.columns if "MB Error" in x][0]
    return float(w[col].iloc[-1])

def qaqc():
    flags = []
    for c in CASES:
        so = open(os.path.join(RUNS, f"{P}_{c}", "raven_stdout.txt")).read()
        qa(f"{c}: Raven reports 'Successful Simulation'", "Exiting Gracefully: Successful Simulation" in so, flags)
    # A1
    r0, r1 = os.path.join(RUNS, f"{P}_D0ref"), os.path.join(RUNS, f"{P}_D0")
    files = sorted(f for f in os.listdir(r0) if f not in ("raven_stdout.txt", "Raven_errors.txt"))
    diff = [f for f in files if not filecmp.cmp(os.path.join(r0, f), os.path.join(r1, f), shallow=False)]
    qa(f"A1 D0 byte-identical to D0ref in {len(files)} output files (differing: {diff})", len(diff) == 0 and len(files) > 0, flags)
    qa("A1 D0 writes no ReachLoss.csv (off by default, R8)", not os.path.exists(os.path.join(r1, f"{STEM}_ReachLoss.csv")), flags)
    # A2
    h0, h1 = rd("D0", "Hydrographs"), rd("D1", "Hydrographs")
    num = [c for c in h0.columns if c.endswith("[m3/s]")]
    mx = float(np.nanmax(np.abs(h0[num].values - h1[num].values)))
    qa(f"A2 D1 Hydrographs vs v2: max |diff| {mx:.3e} m3/s <= 1e-4 over {len(num)} columns, {len(h0)} rows", mx <= 1e-4 and len(h0) == len(h1), flags)
    rl = {c: pd.read_csv(os.path.join(RUNS, f"{P}_{c}", f"{STEM}_ReachLoss.csv")) for c in ["D1", "D2", "D3", "D4", "D5a", "D5b"]}
    d1 = rl["D1"].copy(); d1["date"] = pd.to_datetime(d1["date"])
    wy = np.where(d1.date.dt.month >= 10, d1.date.dt.year + 1, d1.date.dt.year)
    v = d1.assign(wy=wy).query("wy >= 1994 and wy <= 2024").groupby("wy").vol_step_m3.sum() / 1e6   # 1 GL = 1e6 m3
    # E58: F80's 3.271 GL/yr is from bdr_wm_dem (D49) and an offline series, not v2; the D32 law itself is the benchmark
    jd1 = d1.date.dt.dayofyear.values - 1
    lawq = np.where((jd1 >= 166) & (jd1 <= 242), np.minimum(0.17 * d1.Qpre_m3s.values, 0.55), 0.0)
    e = float(np.nanmax(np.abs(lawq - d1.Qact_m3s.values)))
    qa(f"A2 D1 daily loss equals the D32 law min(0.17 Qpre, 0.55) in window, 0 outside (max |diff| {e:.2e} <= 1e-9 m3/s)", e <= 1e-9, flags)
    print(f"A2 info: D1 loss WY1994-2024 {v.mean():.4f} GL/yr ({len(v)} WY); F80 3.271 is on bdr_wm_dem, not v2 (not a test)")
    # A3
    for c in ["D2", "D3"]:
        dmb = mb_end(c) - mb_end("D1")
        qa(f"A3 {c} end-of-run MB error minus D1's: {dmb:.4f} mm (|.| <= 1 mm)", abs(dmb) <= 1.0, flags)
    # A4, A5, A6
    for c, t in rl.items():
        h = rd(c, "Hydrographs")
        qcols = [x for x in h.columns if x.endswith("[m3/s]")]
        neg = int((h[qcols].values < -1e-9).sum())
        qa(f"A4 {c}: no negative flows in Hydrographs ({len(qcols)} columns)", neg == 0, flags)
        q209 = h[[x for x in qcols if x.startswith("Mission_Creek209") and "obs" not in x][0]].values
        # start-of-step (pre-diversion) flow for ReachLoss row k: Hydrographs row k (post-diversion) + actual loss of step k-1
        n = len(t)
        qa(f"A4 {c}: Hydrographs and ReachLoss.csv aligned (rows {len(q209)} vs {n}, first dates {h.date.iloc[0]} vs {t.date.iloc[0]})", len(q209) >= n and h.date.iloc[0] == t.date.iloc[0], flags)
        qpre = t.Qpre_m3s.values   # flow the loss law used (written by the patch, R7)
        # reconciliation, two paths (P1): end-of-step point flow P[k] = Qpre[k+1] - Qact[k]; Hydrographs report step averages
        Pend = qpre[1:] - t.Qact_m3s.values[:-1]
        avg = 0.5 * (Pend[:-1] + Pend[1:])
        best = min(np.nanmax(np.abs(avg - q209[off:off + len(avg)])) for off in (1, 2))
        # tolerance = Hydrographs print precision (6 significant digits): 1e-5 + 5e-6 |Q| (E58)
        best = min(np.nanmax(np.abs(avg - q209[off:off + len(avg)]) - (1e-5 + 5e-6 * np.abs(q209[off:off + len(avg)]))) for off in (1, 2))
        qa(f"A4 {c}: Qpre reconciles with the averaged Hydrographs within print precision (worst excess {best:.2e} <= 0)", best <= 0, flags)
        ok = bool(np.all(t.Qact_m3s.values <= qpre + 1e-9))
        qa(f"A4 {c}: Qact <= start-of-step flow every day ({len(t)} days)", ok, flags)
        vs = 0.5 * (t.Qact_m3s.values[1:] + t.Qact_m3s.values[:-1]) * 86400.0
        rel = np.max(np.abs(vs - t.vol_step_m3.values[1:]) / np.maximum(np.abs(vs), 1.0))
        qa(f"A5 {c}: vol_step equals the trapezoid of Qact (max rel {rel:.2e} <= 1e-6)", rel <= 1e-6, flags)
        if c == "D3":
            s = t.store_m3.values
            pred = s[:-1] * math.exp(-1.0 / 56.0) + t.vol_step_m3.values[1:]   # release at step start, then the new loss
            rel = np.max(np.abs(pred - s[1:]) / np.maximum(s[1:], 1.0))
            qa(f"A6 D3: return store follows the linear-store law k 56 d (max rel {rel:.2e} <= 1e-9)", rel <= 1e-9, flags)
            qa(f"A6 D3: store is never negative (min {s.min():.3f} m3)", s.min() >= 0, flags)
        if c in ("D4", "D5a", "D5b"):
            jd = pd.to_datetime(t.date).dt.dayofyear.values - 1
            inwin = (jd >= 166) & (jd <= 242)
            ref = np.where(inwin, law(c, qpre), 0.0)
            err = np.nanmax(np.abs(ref - t.Qact_m3s.values))
            qa(f"A7 {c}: Raven loss equals the reference law at the start-of-step flow (max |diff| {err:.2e} <= 1e-8 m3/s)", err <= 1e-8, flags)
    a, b = rl["D4"].Qact_m3s.values, rl["D5a"].Qact_m3s.values
    qa("A7 P3 disconnection limit: D5a (water table below bed bottom) equals D4 exactly", np.array_equal(a, b), flags)
    t = rl["D5b"]; gc = int(t.gain_clip_steps.iloc[-1])
    qa(f"A7 P2 no flow against the gradient: D5b gain clipped on {gc} days (> 0) and no loss below stage h_gw", gc > 0, flags)
    # P8 report
    h = rd("D1", "Hydrographs"); h["date"] = pd.to_datetime(h["date"])
    q = h[[x for x in h.columns if x.startswith("Mission_Creek209") and "obs" not in x][0]]
    jd = h.date.dt.dayofyear - 1
    jump = float(np.nanmedian(np.abs(q.diff()[jd == 166]) / q.shift(1)[jd == 166]))
    qa(f"P8 D1 season-boundary step at 16 Jun: median day-to-day change {jump:.1%} (inherited from D32's window)", False, flags, known=True)
    print(f"qaqc: {len(flags)} FLAG" + (": " + "; ".join(flags) if flags else ""))
    sys.exit(1 if flags else 0)

if __name__ == "__main__":
    {"build": build, "run": run, "qaqc": qaqc}[sys.argv[1]]()
