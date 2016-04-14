#!/bin/bash

# extracts number of sites under selection from json output of resultsBusted
# we are looking for evidence rations section, then optimized null
# these values are for each site without gaps
# we return if (2*log2(X) >= 9) print 1; else print 0;
# if site is under selection we print 1 otherwise we print 0;
# finally we print only numbers of sites with positive selection

# for gammarus project

# usage hyphy.json2sites.sh file.json

cat $1  | awk 'BEGIN{flag=0;}{if ($0 ~ /evidence ratios/) flag=1; if (($0 ~ /optimized null/) && (flag==1)) {getline;print $0;}}' | sed s/","/"\n"/g | sed s/"\["// | sed s/"\]"// | awk '{print log($0)/log(2)}' | awk '{if (2*$0>=9) print 1;else print 0;}' | awk '{if ($0 == 1) print NR}'