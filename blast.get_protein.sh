#!/bin/bash

efetch -db protein -id $1 -format fasta | grep .
# | xtract -element "INSDSeq_sequence"