#!/bin/bash

# activate conda env

bname=`basename $3 .cram`
verifybamid2 --SVDPrefix $1 --Reference $2 --BamFile $3 --Output $bname
