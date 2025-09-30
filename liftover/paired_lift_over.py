#!/usr/bin/env python3
# Usage: python paired_liftover.py orig.hg19.vcf lifted.hg38.bed out.hg38.vcf

import sys

if len(sys.argv) != 4:
    print("Usage: paired_liftover_simple.py orig.hg19.vcf lifted.hg38.bed out.hg38.vcf", file=sys.stderr)
    sys.exit(1)

vcf_in, bed_in, vcf_out = sys.argv[1], sys.argv[2], sys.argv[3]

def parse_bed(line):
    c,s,e = line.rstrip("\n").split("\t")[:3]
    return c, int(s), int(e)

def replace_info(info, end, svlen):
    parts = [p for p in info.split(";") if p and not p.startswith("END=") and not p.startswith("SVLEN=")]
    parts.append(f"END={end}")
    parts.append(f"SVLEN={svlen}")
    return ";".join(parts)

with open(vcf_in) as vf, open(bed_in) as bf, open(vcf_out, "w") as out:
    # copy header
    for line in vf:
        if line.startswith("#"):
            out.write(line)
        else:
            # first non-header line found -> rewind file iterator to include this line in records
            first_record = line
            break
    else:
        sys.exit("No VCF records found")

    vlines = [first_record.rstrip("\n")] + [l.rstrip("\n") for l in vf]
    blines = [l.rstrip("\n") for l in bf if l.strip() and not l.startswith("#")]

    n = min(len(vlines), len(blines))
    for i in range(n):
        vcols = vlines[i].split("\t")
        bchrom, bstart, bend = parse_bed(blines[i])
        new_pos = str(bstart + 1)
        # usually for bed len  end-start, but when converting from vcf we should -1
        svlen_abs = bend - bstart - 1
        # preserve sign if original SVLEN negative
        orig_info = vcols[7] if len(vcols) > 7 else "."
        orig_svlen = None
        for t in orig_info.split(";"):
            if t.startswith("SVLEN="):
                try:
                    orig_svlen = int(t.split("=",1)[1]); break
                except: pass
        new_svlen = -abs(svlen_abs) if orig_svlen is not None and orig_svlen < 0 else svlen_abs
        new_info = replace_info(orig_info, bend, new_svlen)
        # ensure there are at least 8 columns
        while len(vcols) < 8:
            vcols.append(".")
        vcols[0] = bchrom
        vcols[1] = new_pos
        vcols[7] = new_info
        out.write("\t".join(vcols) + "\n")
