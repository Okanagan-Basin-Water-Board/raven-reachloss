#!/usr/bin/env bash
# Negative tests for the :ReachLoss parser guards (D181 hardening). Each case puts a deliberately bad block
# into a copy of the D4 short case and must be rejected with the expected message, never run to completion.
# Needs the cases from `reachloss_mc.py build` and the TEST binary. One PASS/FLAG line per case; exits 1 on any FLAG.
set -uo pipefail
MC="${MC:-$HOME/Projects/mission-creek-mass-balance}"
BIN="${BIN:-$(cd "$(dirname "$0")/.." && pwd)/dist/raven-v4.12-reachloss-TEST-macos-arm64/Raven}"
STEM=natural-calibrations-mission-nat-v1-1-0-cal1-validation
SRC="$MC/models/rl181_D4_v412"
W="$(mktemp -d "${TMPDIR:-/tmp}/rl_neg.XXXXXX")"; trap 'rm -rf "$W"' EXIT
fail=0
run_neg() {
  local name=$1 block=$2 expect=$3 d="$W/$1"
  mkdir -p "$d/out"
  for f in "$SRC"/*; do ln -s "$(python3 -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' "$f")" "$d/"; done
  rm "$d/managed_operations.rvt"
  python3 - "$SRC/managed_operations.rvt" "$d/managed_operations.rvt" "$block" <<'EOF'
import sys
t = open(sys.argv[1]).read(); a = t.find(":ReachLoss 209"); e = t.find(":EndReachLoss", a) + len(":EndReachLoss")
assert a >= 0 and e > a, "no :ReachLoss block in the D4 case"
open(sys.argv[2], "w").write(t[:a] + sys.argv[3].replace("\\n", "\n") + t[e:])
EOF
  (cd "$d" && "$BIN" "$STEM" -o out/ > out/stdout.txt 2>&1)
  if grep -q "$expect" "$d/out/stdout.txt" "$d/out/Raven_errors.txt" 2>/dev/null && ! grep -q "Successful Simulation" "$d/out/stdout.txt"; then
    echo "$name: rejected with '$expect' [PASS]"
  else
    echo "$name: not rejected as expected [FLAG]"; fail=1
  fi
}
run_neg wrong_param_count ':ReachLoss 209 DISCONNECTED 166 242\n  :Parameters 0.2 1.0 3000\n  :Destination SINK\n:EndReachLoss' "wrong number of :Parameters"
run_neg exponent_out_of_range ':ReachLoss 209 DISCONNECTED 166 242\n  :Parameters 0.214 1.0 3000.0 8.0 1.5 0.3 0.4 0.0\n  :Destination SINK\n:EndReachLoss' "exponents bW and bH"
run_neg nan_parameter ':ReachLoss 209 EMPIRICAL_Q 166 242\n  :Parameters nan 1.0 0.55\n  :Destination SINK\n:EndReachLoss' "not a number"
run_neg missing_hru_group ':ReachLoss 209 EMPIRICAL_Q 166 242\n  :Parameters 0.17 1.0 0.55\n  :Destination HRUGROUP NoSuchGroup SOIL[2]\n:EndReachLoss' "destination HRU group not found"
run_neg julian_out_of_range ':ReachLoss 209 EMPIRICAL_Q 400 242\n  :Parameters 0.17 1.0 0.55\n  :Destination SINK\n:EndReachLoss' "julian start and end days"
run_neg unknown_form ':ReachLoss 209 LEAKY 166 242\n  :Parameters 0.17 1.0 0.55\n  :Destination SINK\n:EndReachLoss' "form must be EMPIRICAL_Q"
run_neg missing_destination ':ReachLoss 209 EMPIRICAL_Q 166 242\n  :Parameters 0.17 1.0 0.55\n:EndReachLoss' "Destination missing"
[ "$fail" -eq 0 ] && echo "negative_tests: 0 FLAG" || { echo "negative_tests: FLAGGED"; exit 1; }
