#!/bin/bash

cat $1 | egrep "(RYR1|CACNA1S|ASPH)" | grep missense_variant | awk -F "\t" '{if ($37<0.01) print $0}'
