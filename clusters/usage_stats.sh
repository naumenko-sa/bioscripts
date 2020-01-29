#!/bin/bash

# installation
# http://www.leonardoborda.com/blog/how-to-configure-sysstatsar-on-ubuntudebian/

#start=22:21:00
start=$1
end=$2

# CPU load
sar -q -s $start -e $end |  awk '{print $1","$4}'  | sed 1d | sed 1d > cpu.csv

# memory
sar -r -s $start -e $end | awk '{print $5}' | sed 1d | sed 1d > mem.csv

# IO
sar -b -s $start -e $end | awk '{print $5","$6}' | sed 1d | sed 1d > io.csv

paste -d "," cpu.csv mem.csv io.csv > usage.csv

