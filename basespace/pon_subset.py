#!/usr/bin/env python3

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Subset columns in a TSV file based on selected samples")
    parser.add_argument('--pon_combined_counts', required=True, help='Input TSV file with combined counts')
    parser.add_argument('--samples_select', required=True, help='File listing selected sample names')
    parser.add_argument('--pon_selected_counts', required=True, help='Output TSV file with selected columns')
    return parser.parse_args()

def read_sample_list(file):
    with open(file) as f:
        return [line.strip() for line in f.readlines() if line.strip() and not line.startswith("sample_name")]

def main():
    args = parse_args()

    # Read selected samples
    selected_samples = read_sample_list(args.samples_select)

    # Read header lines and detect the actual header
    with open(args.pon_combined_counts) as f:
        lines = f.readlines()

    meta_lines = [line for line in lines if line.startswith("#")]
    data_start_index = len(meta_lines)
    header_line = lines[data_start_index]

    # Read the actual data into a DataFrame
    df = pd.read_csv(args.pon_combined_counts, sep="\t", comment="#", dtype=str)
    #df_columns = df.columns.str.strip().tolist()

    # Columns to keep
    base_columns = ["contig", "start", "stop", "name"]
    selected_columns = [s for s in selected_samples if s in df.columns]

    all_columns = base_columns + selected_columns

    subset_df = df[all_columns]

    # Write the output with the original metadata lines
    with open(args.pon_selected_counts, 'w') as out:
        for line in meta_lines:
            out.write(line)
        subset_df.to_csv(out, sep="\t", index=False)

if __name__ == '__main__':
    main()
