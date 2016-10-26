#!/bin/bash

bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -t ipython -n 400 -s torque -q parallel_long -r walltime=25:00:00 -r vmem=25G --timeout 600 -r minconcores=2
