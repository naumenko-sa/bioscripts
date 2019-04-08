#!/bin/bash

# modeltest.sh file.fasta num_threads

java -jar /hpf/largeprojects/ccmbio/naumenko/tools/jmodeltest-2.1.10/jModelTest.jar -d $1 -s 3 -g 4 -f -BIC -a -tr $2 -o $1.modeltest
#-s 11 -i
