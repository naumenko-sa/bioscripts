#!/bin/bash

qstat -nt1 | grep R | awk '{print $2" " $12}' | grep node | awk -F '+' '{print $1" "NF}' | awk '{print $1" "$3}' | summary_table1.pl