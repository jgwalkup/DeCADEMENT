#!/usr/bin/env bash
#
# delete_partials.sh
# Removes the community directories flagged as "partial" in missing_communities.csv.
# Missing directories are skipped (nothing to delete).
#
# SAFETY: runs as a DRY RUN by default. It prints what it would delete and
# stops. Re-run with --apply to actually remove the directories.
#
# Usage:
#   ./delete_partials.sh            # dry run, lists targets
#   ./delete_partials.sh --apply    # actually delete

set -euo pipefail

GRADIENT="/pub/jwalkup/DecoupledEnergyDEMENT/simulations_20260606/temp_gradient"
CSV="$GRADIENT/missing_communities.csv"

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

[ -f "$CSV" ] || { echo "ERROR: $CSV not found" >&2; exit 1; }

count=0
while IFS=, read -r temp community status nfiles; do
  [ "$status" = "partial" ] || continue
  target="$GRADIENT/$temp/$community"

  if [ ! -d "$target" ]; then
    echo "SKIP (not found): $target"
    continue
  fi

  count=$((count + 1))
  if [ "$APPLY" -eq 1 ]; then
    echo "DELETING: $target  (had $nfiles files)"
    rm -rf "$target"
  else
    echo "WOULD DELETE: $target  (has $nfiles files)"
  fi
done < <(tail -n +2 "$CSV")

echo
if [ "$APPLY" -eq 1 ]; then
  echo "Deleted $count partial director$([ "$count" -eq 1 ] && echo y || echo ies)."
else
  echo "Dry run: $count partial director$([ "$count" -eq 1 ] && echo y || echo ies) would be deleted."
  echo "Re-run with --apply to delete."
fi
