#!/usr/bin/env python3

"""
extract reads from a bam file and a list
write a fasta file
useful benchmark: 
https://timoast.github.io/blog/2015-10-12-extractreads/
"""

import pysam

def extract_reads(options):
    with open(options.names, "r") as f:
        n = f.readlines()
    
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
                    
    f_out = open(options.out, "w")
    for name in n:
        name = name.rstrip()
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                f_out.write(f">{x.query_name}_{x.reference_name}_{x.reference_start+1}_{x.cigarstring}\n")
                f_out.write(x.query_alignment_sequence + "\n")
    
    f_out.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = "Extract reads by read name from the bam (all hits) and write to fasta")
    parser.add_argument("-b", "--bam", help = "bam file", required = True)
    parser.add_argument("-n", "--names", help = "list of read names to extract", required = True)
    parser.add_argument("-o", "--out", help = "output.fasta", required = True)
    options = parser.parse_args()
    extract_reads(options)