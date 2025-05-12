#!/usr/bin/env python

# extracts information about omim inheritance modes from genemap2.txt
# modified from here: https://github.com/OMIM-org/genemap2-parser/blob/master/parseGeneMap2.py

# The OMIM genemap2 file can downloaded from https://omim.org/downloads
# (registration required).

import sys
import re

inheritance_modes = {}

with open("genemap2.txt") as fomim:
    for line in fomim:
        # Skip comments
        if line.startswith('#'):
            continue

    # Strip trailing new line
    line = line.strip('\n')

    # Get the values
    valueList = line.split('\t')

    # Get the fields
    chromosome = valueList[0]
    genomicPositionStart = valueList[1]
    genomicPositionEnd = valueList[2]
    cytoLocation = valueList[3]
    computedCytoLocation = valueList[4]
    mimNumber = valueList[5]
    geneSymbols = valueList[6]
    geneName = valueList[7]
    approvedGeneSymbol = valueList[8]
    entrezGeneID = valueList[9]
    ensemblGeneID = valueList[10]
    comments = valueList[11]
    phenotypeString = valueList[12]
    mouse = valueList[13]

    # Skip empty phenotypes
    if not phenotypeString:
        continue

    if not ensemblGeneID:
        continue

    print(f"{ensemblGeneID}	{approvedGeneSymbol}")

    # Parse the phenotypes
    for phenotype in phenotypeString.split(';'):
        print(f"PHENO: {phenotype}")

        # Clean the phenotype
        phenotype = phenotype.strip()

        # Long phenotype
        matcher = re.match(r'^(.*),\s(\d{6})\s\((\d)\)(|, (.*))$', phenotype)
        if matcher:

            # Get the fields
            phenotype = matcher.group(1)
            phenotypeMimNumber = matcher.group(2)
            phenotypeMappingKey = matcher.group(3)
            inheritances = matcher.group(5)

            # Get the inheritances, may or may not be there
            # molecular basis for the disorder is known, mutations found in the gene
            print(phenotypeMappingKey)
            if phenotypeMappingKey == '3' and inheritances:
                for inheritance in inheritances.split(','):
                    inheritance = inheritance.strip()
                ar_inh = inheritances.split(',')
                for s_i in ar_inh:
                    s_i1 = s_i.strip()
                    if s_i1 in inheritance_modes:
                        inheritance_modes[s_i1]+=1
                    else:
                        inheritance_modes[s_i1]=1

        # Short phenotype
        else:

            # Short phenotype
            matcher = re.match(r'^(.*)\((\d)\)(|, (.*))$', phenotype)
            if matcher:

                # Get the fields
                phenotype = matcher.group(1)
                phenotypeMappingKey = matcher.group(2)
                inheritances = matcher.group(3)

                # Get the inheritances, may or may not be there
                if inheritances:
                    for inheritance in inheritances.split(','):
                        inheritance = inheritance.strip()


# write inheritance dictionary to see the statistics over categories
sorted_inheritance = dict(sorted(inheritance_modes.items(), \
                           key=lambda item: item[1], reverse=True))
with open("inheritance_dictionary.tsv", "w") as f:
    f.write(f"inheritance_mode\tcount\n")
    for key, value in sorted_inheritance.items():
        f.write(f"{key}\t{value}\n")
