#!/usr/bin/env bash
# Build a TEST Raven v4.12 binary with the OBWB reach-loss patch (decision D181), macOS arm64.
# TEST build only: never for production runs or the competition release (method doc Section 6.1).
# The patch is applied to a temporary copy of upstream Raven v4.12 (commit 81011f0); the source tree is never modified.
# Usage:  RAVEN_SRC=<path to RavenHydroFramework checkout at 81011f0> bash scripts/build_reachloss_macos.sh
#         (default RAVEN_SRC=~/Projects/raven-macos/src). Output: dist/raven-v4.12-reachloss-TEST-macos-<arch>/ ; refuses to overwrite.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${RAVEN_SRC:-$HOME/Projects/raven-macos/src}"
PATCH="$ROOT/patches/OBWB_RavenReachLoss_v1_20261031.diff"
PATCH_SHA="75422a53e4640bc4a4847c2d56a8072a5fec52ed7258f442d7b0099cc2a85b9e"
EXPECTED_SHA="81011f061e53ef1289ebcab36dc3622cea05c76d"
ARCH="$(uname -m)"
OUT="${REACHLOSS_OUT:-$ROOT/dist/raven-v4.12-reachloss-TEST-macos-${ARCH}}"
BREW="$(brew --prefix)"
OPT="-O0"   # same profile as the project's reference build, so D-0 can be byte-identical

fail=0
qa() { if eval "$2"; then echo "$1 [PASS]"; else echo "$1 [FLAG]"; fail=1; fi; }
stop() { [ "$fail" -eq 0 ] || { echo "QAQC failed at stage '$1'; stopping."; exit 1; }; }
fsha() { shasum -a 256 "$1" | awk '{print $1}'; }

[ ! -e "$OUT" ] || { echo "$OUT exists; not overwriting"; exit 1; }

echo "== Source checks"
qa "source tree is Raven v4.12 commit $EXPECTED_SHA" "[ \"\$(git -C '$SRC' rev-parse HEAD)\" = '$EXPECTED_SHA' ]"
qa "source tree has no tracked modifications" "[ -z \"\$(git -C '$SRC' status --porcelain --untracked-files=no)\" ]"
qa "patch sha256 is $PATCH_SHA" "[ \"\$(fsha '$PATCH')\" = '$PATCH_SHA' ]"
stop source

TMPD="$(mktemp -d "${TMPDIR:-/tmp}/raven_reachloss.XXXXXX")"
trap 'rm -rf "$TMPD"' EXIT
B="$TMPD/raven"
mkdir -p "$B"
git -C "$SRC" archive HEAD | tar -x -C "$B"
(cd "$B" && patch -s -p1 < "$PATCH")
qa "D181 markers: SubBasin.h 2, SubBasin.cpp 5, Solvers.cpp 2, Model.cpp 2, ParseTimeSeriesFile.cpp 2" \
   "[ \$(grep -c 'OBWB patch D181' '$B/src/SubBasin.h') -eq 2 ] && [ \$(grep -c 'OBWB patch D181' '$B/src/SubBasin.cpp') -eq 5 ] && [ \$(grep -c 'OBWB patch D181' '$B/src/Solvers.cpp') -eq 2 ] && [ \$(grep -c 'OBWB patch D181' '$B/src/Model.cpp') -eq 2 ] && [ \$(grep -c 'OBWB patch D181' '$B/src/ParseTimeSeriesFile.cpp') -eq 2 ]"
qa "only the five patched files differ from upstream" \
   "[ \"\$(diff -rq '$B/src' '$SRC/src' | grep -v '^Only in $SRC' | sed 's|.*/||' | sort | tr '\n' ' ')\" = 'Model.cpp differ ParseTimeSeriesFile.cpp differ Solvers.cpp differ SubBasin.cpp differ SubBasin.h differ ' ]"
stop patch

echo "== Configure ($OPT) and build"
mkdir -p "$B/lib"
ln -s "$BREW/opt/lp_solve/include" "$B/lib/lp_solve_unix"
ln -s "$BREW/opt/lp_solve/lib"     "$B/lib/lp_solve"
cmake -S "$B" -B "$B/build" -DCMAKE_BUILD_TYPE=None -DCMAKE_PREFIX_PATH="$BREW" -DnetCDF_FOUND=ON -DLPSOLVE=ON \
  -DCMAKE_CXX_FLAGS="$OPT -I$BREW/include" -DCMAKE_EXE_LINKER_FLAGS="-L$BREW/lib" >/dev/null
if [ -f "$SRC/build-reference/CMakeFiles/Raven.dir/flags.make" ]; then
  qa "compile flags equal the reference build's flags.make" \
     "diff <(grep -E 'CXX_FLAGS|CXX_DEFINES' '$SRC/build-reference/CMakeFiles/Raven.dir/flags.make') <(grep -E 'CXX_FLAGS|CXX_DEFINES' '$B/build/CMakeFiles/Raven.dir/flags.make') >/dev/null"
fi
cmake --build "$B/build" -j "$(sysctl -n hw.ncpu)" >/dev/null 2>&1
qa "binary contains CSubBasin::AddReachLoss and CalcReachLossRate" "[ \$(nm '$B/build/Raven' | grep -c 'AddReachLoss') -ge 1 ] && [ \$(nm '$B/build/Raven' | grep -c 'CalcReachLossRate') -ge 1 ]"
stop build

echo "== Bundle"
mkdir -p "$OUT/libs"
cp "$B/build/Raven" "$OUT/Raven"
dylibbundler -od -b -x "$OUT/Raven" -d "$OUT/libs" -p @executable_path/libs/ >/dev/null 2>&1
codesign --force -s - "$OUT"/libs/*.dylib "$OUT/Raven" 2>/dev/null
cp "$SRC/LICENSE" "$OUT/LICENSE-Raven.txt"
cp "$PATCH" "$OUT/"
{
  echo "Raven v4.12 for macOS $ARCH, profile 'reference' ($OPT), OBWB reach-loss TEST build (decision D181). NOT FOR PRODUCTION."
  echo "Source: https://github.com/CSHS-CWRA/RavenHydroFramework commit $EXPECTED_SHA MODIFIED by $(basename "$PATCH") sha256 $PATCH_SHA"
  echo "Built: $(date -u +%Y-%m-%dT%H:%M:%SZ), macOS $(sw_vers -productVersion), $(clang --version | head -1)"
  echo "Version string: $("$OUT/Raven" -v 2>&1 | head -1)"
} > "$OUT/BUILD_INFO.txt"

echo "== QAQC"
qa "binary is Mach-O $ARCH" "file '$OUT/Raven' | grep -q '$ARCH'"
qa "no load command points into Homebrew" "! otool -L '$OUT/Raven' '$OUT'/libs/*.dylib | grep -q '$BREW'"
qa "code signatures verify" "codesign -v '$OUT/Raven' && for f in '$OUT'/libs/*.dylib; do codesign -v \"\$f\" || exit 1; done"
qa "version string is 4.12 with netCDF and lp_solve" "'$OUT/Raven' -v 2>&1 | head -1 | grep -Eq '^4.12 w/ netCDF w/ lp_solve'"
qa "source tree still clean after the build" "[ -z \"\$(git -C '$SRC' status --porcelain --untracked-files=no)\" ]"
echo "sha256 reachloss TEST $(fsha "$OUT/Raven")  $OUT/Raven"
[ "$fail" -eq 0 ] && echo "QAQC build_reachloss_macos: 0 FLAG" || { echo "QAQC build_reachloss_macos: FLAGGED"; exit 1; }
