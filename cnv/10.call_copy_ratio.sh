#!/bin/bash

# #1 = T.cr.seg

bname=`basename $1 .cr.seg`

gatk CallCopyRatioSegments \
--input $1 \
--output $bname.called.seg
