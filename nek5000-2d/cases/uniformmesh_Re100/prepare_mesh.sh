#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

NEK_BASE="${NEK_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

find_tool() {
  local name="$1"
  if command -v "$name" >/dev/null 2>&1; then
    command -v "$name"; return 0
  fi
  for p in \
    "$NEK_BASE/bin/$name" \
    "$NEK_BASE/tools/$name" \
    "$NEK_BASE/tools/$name/$name"; do
    if [[ -x "$p" ]]; then echo "$p"; return 0; fi
  done
  return 1
}

GENBOX="$(find_tool genbox)" || { echo "FATAL: genbox not found"; exit 2; }
GENMAP="$(find_tool genmap)" || { echo "FATAL: genmap not found"; exit 2; }

rm -f box.re2 set.re2 set.ma2
printf 'set.box\n' | "$GENBOX"
[[ -f box.re2 ]] || { echo "FATAL: genbox did not create box.re2"; exit 2; }
mv box.re2 set.re2

# Standard genmap workflow: case name followed by vertex tolerance.
printf 'set\n0.01\n' | "$GENMAP"
[[ -f set.ma2 ]] || { echo "FATAL: genmap did not create set.ma2"; exit 2; }

echo "MESH READY: $(pwd)"
