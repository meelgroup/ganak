#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <id>" >&2
  exit 1
fi

id="$1"

if ! [[ "$id" =~ ^[0-9]+$ ]]; then
  echo "error: '$id' is not a number" >&2
  exit 1
fi

find out-ganak-mccomp2324-"${id}"-*/ -type f -print0 \
  | xargs -0 -P 10 -n 1 unxz
