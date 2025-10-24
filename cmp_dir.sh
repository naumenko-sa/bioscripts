#!/usr/bin/env bash
set -euo pipefail

dir1=$1
dir2=$2

# resolve to absolute paths (works with realpath or with a POSIX fallback)
if command -v realpath >/dev/null 2>&1; then
  dir1_abs=$(realpath "$dir1")
  dir2_abs=$(realpath "$dir2")
else
  dir1_abs=$(cd "$dir1" && pwd -P)
  dir2_abs=$(cd "$dir2" && pwd -P)
fi

# ensure no trailing slash
dir1_abs=${dir1_abs%/}
dir2_abs=${dir2_abs%/}

find "$dir1_abs" -type f -print0 | while IFS= read -r -d '' f; do
  # f is an absolute path beginning with dir1_abs
  rel=${f#"$dir1_abs"/}   # relative path inside dir1
  target="$dir2_abs/$rel" # corresponding absolute path in dir2

  #echo "SOURCE: $f"
  #echo "REL   : $rel"

  if [ -f "$target" ]; then
    if ! cmp -s "$f" "$target"; then
      echo "DIFFER: $rel"
      cmp "$f" "$target"
    fi
  else
    echo "MISSING in dir2: $rel"
  fi
done
