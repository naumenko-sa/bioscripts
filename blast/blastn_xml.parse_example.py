#!/usr/bin/env python

import os
import sys
from Bio.Blast import NCBIXML

result_handle = open(sys.argv[1])
blast_records = NCBIXML.parse(result_handle)

families = {}

i=0

for blast_record in blast_records:
    # taking one record one hsp
    try:
        hsp = blast_record.alignments[0].hsps[0]
        i += 1
    except:
        continue
    # hsp.align_length
    # ready query with gaps
    # hsp.query
    # reference with NNN
    # hsp.sbjct
    # UMI position
    m = re.search("N{5,}", hsp.sbjct)
    if not m:
        continue
    start = m.start()
    end = m.end() # end + 1 in fact
    # UMI is aligned with NNNNN
    raw_umi = hsp.query[start:end]
    umi = get_umi(raw_umi, families)
    #umi = raw_umi
    #print(f"HSP: {hsp.sbjct}")
    #print(f"UMI: {umi}")
    seq = hsp.query[0:start]
    ref = hsp.sbjct[0:start]

result_handle.close()
