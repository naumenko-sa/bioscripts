#!/usr/bin/env python

"""
Extracts information about OMIM inheritance modes from a genemap2 file.
modified from here: https://github.com/OMIM-org/genemap2-parser/blob/master/parseGeneMap2.py
The OMIM genemap2 file can downloaded from https://omim.org/downloads
(registration required).
MS Copilot used to optimize the script

4 genes don't have names in genemap2.txt
ENSG00000211592	NA
ENSG00000277734	NA
ENSG00000211893	NA
ENSG00000211899	NA

Usage:
    python script.py [--omim_genemap2 INPUT] [--omim_inheritance OUTPUT] [--omim_inheritance_dictionary DICT_OUTPUT]

Defaults:
    --omim_genemap2:                genemap2.txt
    --omim_inheritance:             omim_inheritance.tsv
    --omim_inheritance_dictionary:  omim_inheritance_dictionary.tsv

Output:
   omim_inheritance.tsv: ensembl_gene_id	gene_name	inheritance modes(;)
   omim_inheritance_dictionary: inheritance_mode	count

"""

import re
import os
import csv
import argparse

def update_inheritance(inheritance_str, counts):
    """ Updates the counts dictionary given a comma-separated string of inheritances. """
    for inh in [item.strip() for item in inheritance_str.split(',')]:
        if inh:
            counts[inh] = counts.get(inh, 0) + 1


def process_genemap2(input_file):
    """
    Process the OMIM genemap2 file and return:
      - a list of gene inheritance records (ensembl_geneid, gene_name, omim_inheritance)
      - a dictionary of inheritance mode statistics.
    """
    # Precompile regex patterns for efficiency.
    # The long phenotype pattern captures:
    #   group(1): phenotype description
    #   group(2): 6-digit mim number
    #   group(3): mapping key
    #   group(4): optional inheritances (after a comma and space)
    #   there could be multiple modes in one line: mode1, mode2
    long_pattern = re.compile(
        r'^(.*),\s(\d{6})\s\((\d)\)(?:, (.*))?$'
    )
    # The short phenotype pattern:
    #   group(1): phenotype description
    #   group(2): mapping key
    #   group(3): optional inheritances
    #   Short pattern does not have inheritance info in fact
    #   SHORT: Hemolytic anemia due to phosphofructokinase deficiency (1)
    short_pattern = re.compile(
        r'^(.*)\((\d)\)(?:, (.*))?$'
    )

    inheritance_modes = {}
    # as of 2025, OMIM has a duplicate record for 
    # NSG00000185960	SHOX	Pseudoautosomal dominant;Pseudoautosomal recessive
    genes_present = set()  
    inheritance_records = []

    # Define field names in the order they appear in genemap2.txt.
    fieldnames = [
        "chromosome", "genomic_position_start", "genomic_position_end",
        "cyto_location", "computed_cyto_location", "mim_number", "gene_symbols",
        "gene_name", "approved_gene_symbol", "entrez_geneid", "ensembl_geneid",
        "comments", "phenotype_string", "mouse"
    ]

    with open(input_file, newline='') as f:
        # Filter out any comment lines (lines starting with '#')
        filtered_lines = (line for line in f if not line.startswith('#'))
        reader = csv.DictReader(filtered_lines, delimiter="\t", fieldnames=fieldnames)

        for record in reader:
            phenotype_field = record.get("phenotype_string", "").strip()
            ensembl_geneid = record.get("ensembl_geneid", "").strip()
            gene_name = record.get("approved_gene_symbol", "").strip() or "NA"

            if not phenotype_field or not ensembl_geneid:
                continue

            if ensembl_geneid in genes_present:
                continue
            genes_present.add(ensembl_geneid)

            inheritance4gene = []

            # Each record may contain multiple phenotypes separated by ';'
            for phenotype in map(str.strip, phenotype_field.split(';')):
                match = long_pattern.match(phenotype)
                if match:
                    # Only consider phenotypes with mapping key '3' and provided inheritance info.
                    # (3): molecular basis is known, mutations found in the gene
                    if match.group(3) == '3' and match.group(4):
                        inh_str = match.group(4).strip()
                        update_inheritance(inh_str, inheritance_modes)
                        for inh in map(str.strip, inh_str.split(',')):
                            if inh and inh not in inheritance4gene:
                                inheritance4gene.append(inh)
                else:
                    # Short phenotype pattern matching is ignored per original logic.
                    if short_pattern.match(phenotype):
                        continue

            s_inheritance = ';'.join(sorted(inheritance4gene)) if inheritance4gene else "NA"
            inheritance_records.append((ensembl_geneid, gene_name, s_inheritance))
    
    return inheritance_records, inheritance_modes

def write_inheritance_output(filename, records):
    """Write the gene inheritance output with header."""
    with open(filename, "w") as outf:
        writer = csv.writer(outf, delimiter="\t", lineterminator="\n")
        # Write header: ensembl_gene_id, gene_name, omim_inheritance
        writer.writerow(["ensembl_gene_id", "gene_name", "omim_inheritance"])
        writer.writerows(records)

def write_inheritance_dictionary(filename, modes):
    """Write the inheritance dictionary (stats) to a TSV file."""
    sorted_inheritance = sorted(modes.items(), key=lambda item: item[1], reverse=True)
    with open(filename, "w") as outf:
        writer = csv.writer(outf, delimiter="\t", lineterminator="\n")
        writer.writerow(["inheritance_mode", "count"])
        for mode, count in sorted_inheritance:
            writer.writerow([mode, count])

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract OMIM inheritance info from a genemap2 file."
    )
    parser.add_argument(
        "--omim_genemap2",
        default="genemap2.txt",
        help="Input OMIM genemap2 file (default: genemap2.txt)"
    )
    parser.add_argument(
        "--omim_inheritance",
        default="omim_inheritance.tsv",
        help="Output gene inheritance file (default: omim_inheritance.tsv)"
    )
    parser.add_argument(
        "--omim_inheritance_dictionary",
        default="omim_inheritance_dictionary.tsv",
        help="Output inheritance dictionary file (default: omim_inheritance_dictionary.tsv)"
    )
    args = parser.parse_args()
    return parser, args

def main():
    parser, args = parse_arguments()
    if not os.path.exists(args.omim_genemap2):
        parser.error(f"Input file '{args.omim_genemap2}' not found.")
    records, modes = process_genemap2(args.omim_genemap2)
    write_inheritance_output(args.omim_inheritance, records)
    write_inheritance_dictionary(args.omim_inheritance_dictionary, modes)

if __name__ == "__main__":
    main()
