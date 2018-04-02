#!/usr/bin/env python2
import sys
import collections

lines = sys.stdin.read().splitlines()
counter = collections.Counter(lines)

for key,count in counter.most_common():
      print(str(key)+'\t'+str(1.0*count/len(lines)))
