#!/bin/bash

samtools stats $bam | grep "^SN" | cut -f 2- | head -n1 | awk '{print $4}' > ${bam}.reads_number
