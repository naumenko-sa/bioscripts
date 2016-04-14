#!/bin/bash

samtools stats $1 | grep "^SN" | cut -f 2- | head -n1 | awk '{print $4}'
