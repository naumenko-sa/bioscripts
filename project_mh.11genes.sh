#!/bin/bash

#extracts variants for 11 genes of RYR1 interactome from gzipped vcf file

gunzip -c $1 | egrep "(RYR1|CACNA1S|ASPH|STAC3|TRDN|JPH2|CASQ1|ATP2A1|ATP2A2|CALM1|FKBP1A)" | grep missense_variant
