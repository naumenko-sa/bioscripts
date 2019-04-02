#!/bin/bash

raxml=/hpf/largeprojects/ccmbio/naumenko/tools/standard-RAxML/raxmlHPC-PTHREADS-SSE3
$raxml -T 5 -s $1 -n $1.tree -m GTRGAMMAI -p 8593771075132634533
#-m PROTGAMMAWAG