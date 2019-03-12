#!/bin/bash

/home/tools/RAxML-8.0.0/raxmlHPC-PTHREADS-SSE3 -T 5 -s $1 -n $1.tree -m GTRGAMMAI -p 8593771075132634533
#-m PROTGAMMAWAG