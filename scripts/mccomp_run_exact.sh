#/usr/bin/env bash
file=$1
if [ -z "$file" ]; then
    echo "Usage: $0 <path_to_mccomp_run_file>"
    exit 1
fi

if grep -qE '^c t (mc|pmc)' "$file"; then
    echo "c o The file seems to be unweighted. Will run with --mode 0."
    ./ganak --mode 0 --maxcache=22000 "$file"
elif grep -qE '^c t (wmc|pwmc)' "$file"; then
    echo "c o The file seems to be weighted. Will run with --mode 1."
    ./ganak --mode 1 --maxcache=16000 "$file"
elif grep -qE '^c t (cpxmc|pcpxmc)' "$file"; then
    echo "c o The file seems to be complex number-weighted. Will run with --mode 2."
    ./ganak --mode 2 --maxcache=10000 "$file"
else
  echo "c o The file does not contain a valid mccomp header."
  exit 1
fi
