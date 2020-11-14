#!/usr/bin/python

# merges 2 fastq files for left and right read into single file
# usage: fastq2to1.py left.fq right.fq single.fq

import sys

f1 = open(sys.argv[1])
f2 = open(sys.argv[2])
fout = open(sys.argv[3],'w')

f1_buf = []
f2_buf = []

line = f1.readline()
while(line != ''):
    f1_buf.append(line)
    line = f2.readline()
    f2_buf.append(line)
    for i in range(3):
	f1_buf.append(f1.readline())
	f2_buf.append(f2.readline())
	
    for i in range(4):
	fout.write(f1_buf[i])
    for i in range(4):
	fout.write(f2_buf[i])

    del f1_buf[:]
    del f2_buf[:]
    
    line = f1.readline()

f1.close()
f2.close()
fout.close()
