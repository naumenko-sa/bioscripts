#!/bin/bash

bcbio_nextgen.py ../config/c4r_multisample.yaml -t ipython -n 800 -s torque -q parallel_long -r walltime=200:00:00 -r vmem=50G --timeout 600 -r minconcores=2
