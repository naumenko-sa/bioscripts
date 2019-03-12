#!/bin/bash

#checks conservative 1st and 3rd positions in the codon + 2nd nonconvervative

cat $1 | egrep "^[0-9]{1,3}\s*([ATGC])\1{7}" | egrep "([ATGC])\1{7}$" | awk '{print $1"\t"$3}' | egrep -v "([ATGC])\1{7}" > `echo $1 | sed s/nondeg/2nd_noncons/`;
