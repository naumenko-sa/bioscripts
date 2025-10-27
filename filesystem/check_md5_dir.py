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
from pathlib import Path
import argparse
import subprocess
import sys

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

def verify_md5s(root: Path):
    cmd = "md5sum"

    md_files = sorted(root.rglob("*.md5sum"))
    if not md_files:
        print("No .md5sum files found.")
        return 0

    for mdfile in md_files:
        try:
            expected = read_expected_checksum(mdfile)
            if not expected:
                print(f"ERROR  : {mdfile} : no checksum found in file")
                continue

            target = mdfile.with_name(mdfile.stem)  # sibling file without .md5sum
            if not target.exists():
                print(f"MISSING: {target} (side-by-side with {mdfile.name})")
                continue

            actual = compute_md5_with_cmd(target)
            if actual == expected:
                print(f"OK     : {target}")
            else:
                print(f"MISMATCH: {target} (expected {expected}, got {actual})")
        except Exception as e:
            print(f"ERROR  : {mdfile} : {e}")

    return 0

def main():
    p = argparse.ArgumentParser(description="Verify sibling MD5 checksums using system md5 tool.")
    p.add_argument("root", nargs="?", default=".", help="Root directory to search (default: current dir)")
    args = p.parse_args()
    rc = verify_md5s(Path(args.root))
    sys.exit(rc)

if __name__ == "__main__":
    main()
