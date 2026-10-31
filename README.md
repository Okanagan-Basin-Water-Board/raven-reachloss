# raven-reachloss

A TEST patch for [Raven](https://github.com/CSHS-CWRA/RavenHydroFramework) v4.12 that adds a native **losing-reach (channel transmission loss) process**, with a Mission Creek (Okanagan, BC) test harness. Developed by the Okanagan Basin Water Board (OBWB), Water Science and Stewardship.

> **Status: TEST build.** Not for production runs. The patch is a working prototype for discussion with the Raven developers. It is not an upstream feature.

## Why

Raven v4.12 has no loss term in a river reach (`RiverReach.cpp`). Projects represent a losing reach with `:FlowDiversionLookupTable SBID -1`, which sends the water out of the model. On the Mission Creek alluvial fan, the groundwater model shows two regimes: the reaches above the WSC gauge 08NM116 are disconnected from the aquifer, and the reaches below it are connected. The lost water recharges the fan aquifer rather than leaving the basin. This patch represents both regimes and keeps the water in the watershed balance.

## What the patch adds

A new block in the `.rvt` (time series) file:

```
:ReachLoss [SBID] [EMPIRICAL_Q | DISCONNECTED | CONNECTED] {start_julian_day} {end_julian_day}
  :Parameters [p1] [p2] ...
  :Destination SINK
  :Destination HRUGROUP [HRU group name] [state variable, e.g. SOIL[2]]
  :Destination SUBBASIN [SBID] [k, days]
:EndReachLoss
```

| Form | `:Parameters` | Loss rate L [m³/s], Q = reach outflow at the start of the step |
|---|---|---|
| `EMPIRICAL_Q` | a, b, cap | min(a Q^b, cap) |
| `DISCONNECTED` | Kbed [m/d], d [m], length [m], aW, bW, aH, bH, cap | Kbed (h + d)/d W length, with W = aW Q^bW and h = aH Q^bH |
| `CONNECTED` | as DISCONNECTED, then z_bed [m], h_gw [m] | Kbed ((z_bed + h) - max(h_gw, z_bed - d))/d W length |

Rules:
- cap ≤ 0 means no cap. The loss is never negative and never more than Q (supply limit).
- `CONNECTED` falls back to the disconnected limit as soon as h_gw drops below the bed bottom, so there is no jump at the switch. When h_gw is above stage, the would-be gain is clipped to 0 and counted (`gain_clip_steps`). Gaining reaches need a coupled groundwater store, which is future work.

Destinations:
- `SINK`: the water leaves the model, as `:FlowDiversionLookupTable ... -1` does today.
- `HRUGROUP`: the water is added to a state variable of the group's enabled HRUs, area-weighted. It stays in the watershed balance.
- `SUBBASIN`: the water goes into a linear store with time constant k (days), which releases into the target subbasin as lateral inflow. The store is counted in channel storage.

Each day, `{RunName}_ReachLoss.csv` in the output folder (named by Raven's `FilenamePrepare`, like every other output) records, per reach: the start-of-step flow the loss law used (`Qpre_m3s`), the actual loss (`Qact_m3s`), the step volume, the return-store volume and the gain-clip count. Raven's Hydrographs.csv reports step-averaged flows, so `Qpre_m3s` is the only direct record of the flow the law saw. A model without `:ReachLoss` runs exactly as unpatched Raven.

Input checks and run-time guards:
- The parser stops on: an unknown form or destination; the wrong parameter count; a NaN parameter; a negative Kbed, aW, aH, a or b; bed thickness or length ≤ 0; hydraulic-geometry exponents outside 0 to 1; Julian days outside 0 to 366; a missing HRU group, state variable or subbasin; a destination that is not a water store; a missing `:EndReachLoss`.
- The Julian window wraps across 31 Dec when start > end, as for `:FlowDiversion`.
- At run time the patch stops on: a NaN loss rate; an HRU group with no enabled area; an output file that cannot be opened.
- It warns when constituent transport is active (lost water carries no mass or heat) and under management optimization (the loss is outside the optimization problem).

Changed files: `SubBasin.h`, `SubBasin.cpp`, `Solvers.cpp`, `Model.cpp`, `ParseTimeSeriesFile.cpp`, about 300 lines, each marked `// OBWB patch D181 (reach loss)`.

## Continuous integration

`.github/workflows/build.yml` runs on every push. It:
- checks the patch hash against the build script;
- applies the patch to the pinned upstream commit and checks the marker counts;
- regenerates the edits from `tools/apply_reachloss.py` and diffs them against the patched source;
- compiles on Ubuntu 24.04 with netCDF.

## Build (macOS arm64)

```bash
git clone https://github.com/CSHS-CWRA/RavenHydroFramework raven-src && git -C raven-src checkout 81011f061e53ef1289ebcab36dc3622cea05c76d
RAVEN_SRC=$PWD/raven-src bash scripts/build_reachloss_macos.sh
```

The script applies the patch to a temporary copy, never edits your checkout, prints one PASS/FLAG line per check and stops on any FLAG. It needs Homebrew `netcdf`, `lp_solve`, `cmake` and `dylibbundler`. Output: `dist/raven-v4.12-reachloss-TEST-macos-arm64/`.

`tools/apply_reachloss.py` regenerates the same edits from source, and `patches/` holds the resulting diff.

## Mission Creek test harness

`mission_creek/reachloss_mc.py build | run | qaqc` builds eight cases from the OBWB Mission Creek working model (v2). This needs the private project repository and its forcing data. The cases are:
- D0 and D0ref: v2 on the patched and reference binaries, checking byte-identity;
- D1: D32 rewritten as `:ReachLoss`;
- D2: the loss sent to the fan HRUs' deep soil;
- D3: a 56-day return store to the lake-outlet subbasin;
- D4, D5a and D5b: the physical forms.

The acceptance rules (A1 to A7) are in the script header, written before any run. They test mass conservation, the supply limit, the linear-store law, the loss law, the disconnection limit and the gradient sign.

## Licence

The patch modifies Raven, which is distributed under the Artistic License 2.0 (`LICENSE`). Contributions upstream are subject to the Raven contributor licensing agreement.

## Citation and contact

Okanagan Basin Water Board (2026). raven-reachloss: a native losing-reach process for Raven (TEST). Contact: Nelson Jatel, Senior Director, Water Science and Stewardship, OBWB.

Method: OBWB Mission Creek project, *Losing reaches on the Mission Creek alluvial fan: method and approach for phases A to D* (v1, 2026-10-31), decision D181.
