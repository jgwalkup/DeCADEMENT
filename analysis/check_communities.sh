#!/usr/bin/env bash
#
# check_communities.sh
# Identifies communities that are missing or only partially saved across all
# temperature directories in a temperature-gradient simulation run.
#
# Output: missing_communities.csv with columns: temperature,community,status,n_files
#   status = "missing"  -> community directory does not exist
#   status = "partial"  -> directory exists but has fewer files than a complete run
#
# Usage:
#   ./check_communities.sh [GRADIENT_DIR]
# If GRADIENT_DIR is omitted, the default below is used.

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG -- adjust these to match your run
# ---------------------------------------------------------------------------
GRADIENT="${1:-/pub/jwalkup/DecoupledEnergyDEMENT/DecoupledEnergyDEMENTUpload}"
TEMP_GLOB="static_T*_Pn0p5"   # temperature subdirectory pattern
FIRST=1                       # first community number
LAST=1000                     # last community number
SEED="s42"                    # seed suffix
PAD=4                         # zero-padding width (set to 0 for no padding)
OUT="$GRADIENT/missing_communities.csv"

# ---------------------------------------------------------------------------
# Build expected community-name list
# ---------------------------------------------------------------------------
EXPECTED="$(mktemp)"
if [ "$PAD" -gt 0 ]; then
  seq "$FIRST" "$LAST" \
    | xargs -I{} printf "community_%0${PAD}d_${SEED}\n" {} \
    | sort > "$EXPECTED"
else
  seq "$FIRST" "$LAST" \
    | xargs -I{} printf "community_%d_${SEED}\n" {} \
    | sort -V > "$EXPECTED"
fi

# ---------------------------------------------------------------------------
# Determine the "complete" file count (the mode across all present dirs).
# Files are counted recursively so nested output (per-pulse dirs, Zarr stores)
# is handled correctly.
# ---------------------------------------------------------------------------
COUNTS="$(mktemp)"   # lines: "<nfiles> <tempdir> <community>"
for tdir in "$GRADIENT"/$TEMP_GLOB; do
  [ -d "$tdir" ] || continue
  tname="$(basename "$tdir")"
  for cdir in "$tdir"/community_*_"$SEED"/; do
    [ -d "$cdir" ] || continue
    n="$(find "$cdir" -type f | wc -l)"
    printf "%s %s %s\n" "$n" "$tname" "$(basename "$cdir")" >> "$COUNTS"
  done
done

if [ ! -s "$COUNTS" ]; then
  echo "ERROR: no community directories matched under $GRADIENT/$TEMP_GLOB" >&2
  echo "Check GRADIENT, TEMP_GLOB, and SEED." >&2
  exit 1
fi

COMPLETE="$(awk '{print $1}' "$COUNTS" | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')"
echo "Detected complete file count (mode): $COMPLETE"
echo "File-count distribution (count_of_dirs  n_files):"
awk '{print $1}' "$COUNTS" | sort -n | uniq -c

# ---------------------------------------------------------------------------
# Write combined report: missing dirs + partial dirs
# ---------------------------------------------------------------------------
echo "temperature,community,status,n_files" > "$OUT"

for tdir in "$GRADIENT"/$TEMP_GLOB; do
  [ -d "$tdir" ] || continue
  tname="$(basename "$tdir")"

  # present community names in this temperature dir
  ls -1 "$tdir" 2>/dev/null | sort > "$EXPECTED.present"

  # MISSING: expected but not present
  comm -23 "$EXPECTED" "$EXPECTED.present" | while read -r c; do
    echo "${tname},${c},missing,0" >> "$OUT"
  done

  # PARTIAL: present but fewer files than complete
  awk -v t="$tname" -v c="$COMPLETE" \
      '$2==t && $1 < c {print t","$3",partial,"$1}' "$COUNTS" >> "$OUT"
done

rm -f "$EXPECTED" "$EXPECTED.present" "$COUNTS"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
echo "Wrote: $OUT"
echo "Total flagged: $(($(wc -l < "$OUT") - 1))"
echo "By temperature and status:"
tail -n +2 "$OUT" | cut -d, -f1,3 | sort | uniq -c
