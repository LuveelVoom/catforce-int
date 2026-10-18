#!/bin/bash
# usage: suite.sh <binary> <label> ; per-test clean dirs, dedup off + solution list for exact comparison
# Compare the printed md5s with tests/base.txt. Results go to $OUT_DIR (default: /tmp).
BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); S=$(cd "$(dirname "$0")" && pwd); OUT=${OUT_DIR:-/tmp}/catforce_res_$2; rm -rf "$OUT"; mkdir -p "$OUT"
now(){ perl -MTime::HiRes=time -e 'printf "%.3f", time'; }
if command -v timeout >/dev/null; then TO="timeout 280"; elif command -v gtimeout >/dev/null; then TO="gtimeout 280"; else TO=""; fi
for f in intb conway baitA quad sym ex_orfilter ex_match ex_phase-changing-blinker ex_p34 layout; do
  d=$OUT/$f; mkdir -p "$d"; cp "$S/$f.in" "$d/"; [ -f "$S/$f.rle" ] && cp "$S/$f.rle" "$d/"; cd "$d"; printf "keep-duplicates yes\noutput-list list.txt\n" >> $f.in
  s=$(now); if [ -n "$TO" ]; then $TO "$BIN" $f.in > log.txt 2>&1; else "$BIN" $f.in --max-time 280 > log.txt 2>&1; fi; e=$(now)
  printf "%-26s %6.1fs %-30s %s\n" $f $(echo "$e - $s" | bc) "$(grep '^results' log.txt | tail -1 | cut -d, -f1)" "$(grep -v '^#' list.txt 2>/dev/null | sed 's/^[0-9]*: //' | sort | md5sum | cut -c1-10)"
done
