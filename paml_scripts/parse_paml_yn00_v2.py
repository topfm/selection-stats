#!/usr/bin/env python3

import argparse
import os
import csv
from Bio.Phylo.PAML import yn00

################################################################################
# This script parses PAML yn00 output and outputs a txt file for visualization in R
################################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = f"{filename} is not a file"
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Parse PAML yn00 output')
    parser.add_argument("paml", help="yn00_output file",
                        action=FullPaths,
                        type=is_file)
    return parser.parse_args()

def get_data(infile):
    """Read yn00 output file"""
    return yn00.read(infile)

def write_results(results, out_file):
    """Flatten and write yn00 results to tab-delimited text file"""
    rows = []
    
    for iso in results:
        for comp_iso in results[iso]:
            row = {
                "isolate": iso,
                "comparedTo_iso": comp_iso
            }
            for method in results[iso][comp_iso]:
                for stat_key, stat_val in results[iso][comp_iso][method].items():
                    col_name = f"{method}_{stat_key}"
                    row[col_name] = stat_val
            rows.append(row)

    # Collect all unique headers across all rows
    all_headers = sorted(set().union(*(row.keys() for row in rows)))

    # Write to output file
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

# Main script execution
if __name__ == "__main__":
    args = get_args()
    results = get_data(args.paml)
    
    # Generate output filename
    outWrite = args.paml.rsplit(".", 1)[0] + ".parsed_yn00.txt"
    
    write_results(results, outWrite)
