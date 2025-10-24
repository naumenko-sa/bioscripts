#!/usr/bin/env bash

set -euo pipefail

dir1=$1
dir2=$2

# normalize dir1 to remove any trailing slash for consistent prefix removal
dir1=${dir1%/}
# use find with -print0 and a null-safe read loop to handle spaces/newlines in names
find "./$dir1" -type f -print0 | while IFS= read -r -d '' f; do
  # f is like ./dir1/path/to/file or ./dir1/file
  # remove leading ./ if present
  f_nodot=${f#./}
  # remove the "dir1/" prefix to get the relative path
  rel=${f_nodot#"$dir1"/}

  #echo "SOURCE: $f"
  #echo "REL   : $rel"

  target="${dir2%/}/$rel"  # ensure dir2 has no trailing slash
  if [ -f "$target" ]; then
    if ! cmp -s "$f" "$target"; then
      echo "DIFFER: $rel"
      cmp "$f" "$target"
    fi
  else
    echo "MISSING in dir2: $rel"
  fi
done