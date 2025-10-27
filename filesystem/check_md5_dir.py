#!/usr/bin/env python3
"""
check_md5_dir.py

Recursively find files named *.md5sum. For each one, expect a sibling file
with the same name but without the .md5sum suffix. Read the expected checksum
(from the first token of the first non-empty line) and compute/verify the MD5
by invoking the system md5 tool (md5sum or md5). Prints OK, MISSING, MISMATCH,
or ERROR for each file and continues processing all files.
Exit code 0 on normal completion (does not fail on mismatches).
"""

import sys
import argparse
import subprocess
from pathlib import Path

def read_expected_checksum(mdfile: Path):
    for line in mdfile.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        return line.split()[0].lower()
    return None

def compute_md5_with_cmd(filepath: Path):
    cmd = "md5sum"
    # md5sum outputs: "<hex>  filename"
    proc = subprocess.run([cmd, str(filepath)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"{cmd} failed with code {proc.returncode}")
    return proc.stdout.split()[0].lower()    


def iter_md5_files_with_find(root):
    proc = subprocess.Popen(
        ['find', str(root), '-type', 'f', '-name', '*.md5sum', '-print0'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    assert proc.stdout is not None

    buf = b''
    try:
        while True:
            chunk = proc.stdout.read(8192)
            if not chunk:
                break
            buf += chunk
            while True:
                try:
                    nul_index = buf.index(b'\0')
                except ValueError:
                    break
                raw = buf[:nul_index]
                buf = buf[nul_index+1:]
                if raw:
                    yield Path(raw.decode('utf-8', 'surrogateescape'))
        # any trailing data without NUL (shouldn't happen with -print0) handled here
        if buf:
            yield Path(buf.decode('utf-8', 'surrogateescape'))
    finally:
        proc.stdout.close()
        proc.wait()

def verify_md5s(root: Path):    

    for mdfile in iter_md5_files_with_find(root):
        try:
            expected = read_expected_checksum(mdfile)
            if not expected:
                print(f"ERROR  : {mdfile} : no checksum found in file", flush = True)
                continue

            target = mdfile.with_name(mdfile.stem)  # sibling file without .md5sum
            if not target.exists():
                print(f"MISSING: {target} (side-by-side with {mdfile.name})", flush = True)
                continue

            actual = compute_md5_with_cmd(target)
            if actual == expected:
                print(f"OK     : {target}", flush = True)
            else:
                print(f"MISMATCH: {target} (expected {expected}, got {actual})", flush = True)
        except Exception as e:
            print(f"ERROR  : {mdfile} : {e}", flush = True)           

    return 0

def main():
    p = argparse.ArgumentParser(description="Verify sibling MD5 checksums using system md5 tool.")
    p.add_argument("root", nargs="?", default=".", help="Root directory to search (default: current dir)")
    args = p.parse_args()
    rc = verify_md5s(Path(args.root))
    sys.exit(rc)

if __name__ == "__main__":
    main()
