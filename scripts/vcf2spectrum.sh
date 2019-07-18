#/bin/bash

bcftools stats $1 > $1.stats
cat $1.stats | grep "^ST"  | grep -v "^#" | awk '{print $3","$4}' > $1.spectrum.csv
