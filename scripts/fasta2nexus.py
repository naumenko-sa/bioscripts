#!/usr/bin/env python
import os
import sys

 
if ("--help" in sys.argv) or ("-?" in sys.argv):
    sys.stderr.write("usage: fasta-to-nexus.py [<fasta-file-path>] [<nexus-file-path>]\n")
    sys.exit(1)
 
if len(sys.argv) < 2:
    src = sys.stdin
else:
    src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
    if not os.path.exists(src_fpath):
        sys.stderr.write('Not found: "%s"' % src_fpath)
    src = open(src_fpath)        
 
if len(sys.argv) < 3:
    dest = sys.stdout # os.path.splitext(src_fpath)[0] + ".nex"
else:
    dest_fpath = os.path.expanduser(os.path.expandvars(sys.argv[2]))
    dest = open(dest_fpath, "w")


def writeln(*lines):
    dest.write("\n".join(lines))
    dest.write("\n")


seqs = {}
cur_seq = None
lines = src.readlines()
for i in lines:
    i = i.rstrip("\n\r")  # Strip "\n" and/or "\r" characters from the end.
    if i:
        if i.startswith(">"):
            label = i[1:]
            cur_seq = []
            seqs[label] = cur_seq
        else:
            if cur_seq is None:
                raise Exception("Sequence data found before label")
            cur_seq.extend(i.replace(" ", ""))
 
taxlabels = seqs.keys()
taxlabels.sort()
 
writeln("#NEXUS",
        "",
        "Begin Taxa;",
        "  dimensions ntax=%d;" % len(seqs),
        "  taxlabels")
for taxlabel in taxlabels:
    writeln("    %s" % taxlabel)
writeln("  ;",
        "End;",
        "")
nchar = max(len(s) for s in seqs.values())
writeln("Begin Characters;",
        "  dimensions nchar=%d;" % nchar,
        "  format datatype=dna missing=? gap=-;",
        "  matrix")
for taxlabel in taxlabels:
    writeln("    %s      %s" % (taxlabel, "".join(seqs[taxlabel])))
writeln("  ;",
        "end;")
