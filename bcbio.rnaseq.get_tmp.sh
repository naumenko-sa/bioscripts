#!/bin/bash

ens=`cat annotated_combined.counts | awk -v gene=$1 '{if($3 == gene)print $1}'`
echo $1 `cat combined.gene.sf.tpm | grep $ens`
