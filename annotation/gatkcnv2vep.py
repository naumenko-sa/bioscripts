#!/usr/bin/env python3

# convert gatcnv called.igv.seg output to a format suited for VEP annotation
# http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#sv

import sys
with open(sys.argv[1]) as gatk_seg:
    for line in gatk_seg:
        if not line.startswith("Sample"):
            fields = line.strip().split()
            segment_mean = abs(float(fields[6]))
            if segment_mean >= 2:
                if fields[5] == "+":
                    sv_type = "DUP"
                else:
                    sv_type = "DEL"
                print("\t".join([fields[1], fields[2], fields[3], sv_type, fields[5], fields[4], fields[6]]))

