#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import PAML
from Bio.Phylo.PAML import yn00

################################################################################
# This script parses yn00 output and outputs txt file for visualization in R
################################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='parse paml yn00 output')
    parser.add_argument("paml", help="yn00_output",
                        action=FullPaths,
                        type=is_file)
    return parser.parse_args()

def get_data(infile):
    results = yn00.read(infile)
    return results

def write_results(results):
    with open(outWrite, "w") as outfile:
        outfile.write("isolate\tcomparedTo_iso\tmethod_1\tsanity_check\tYN00_dN\t\tYN00_dN_SE\t\tYN00_kappa\t\tYN00_dS_SE\t\tYN00_N\t\tYN00_S\t\tYN00_t\t\tYN00_omega\t\tYN00_dS\tmethod_2\t\tLWL85m_dN\t\tLWL85m_w\t\tLWL85m_N\t\tLWL85m_S\t\tLWL85m_rho\t\tLWL85m_dS\tmethod3\t\tLPB93_dN\t\tLPB93_dS\t\tLPB93_w\tmethod4\t\tLWL85_dN\t\tLWL85_S\t\tLWL85_dS\t\tLWL85_w\t\tLWL85_N\tmethod5\t\tNG86_dN\t\tNG86_omega\t\tNG86_dS\n")
        for iso in results:
            print(iso)
            for comparison_iso in results[iso]:
                outfile.write("\n"+iso+"\t"+comparison_iso)
                print("compared to "+comparison_iso)
                for stat in results[iso][comparison_iso]:
                    stat_method = stat
                    outfile.write("\t"+stat_method)
                    print(stat)
                    for nested_stat in results[iso][comparison_iso][stat]:
                        method_nested_stat = nested_stat
                        method_nested_stat_value = str(results[iso][comparison_iso][stat][nested_stat])
                        outfile.write("\t"+stat_method+"_"+method_nested_stat+"\t"+method_nested_stat_value)
                        print("nested stat "+nested_stat+" = "+str(results[iso][comparison_iso][stat][nested_stat]))
    return        

args = get_args()
results = get_data(args.paml)
outWrite = args.paml.replace(".output.txt", ".parsed_yn00.txt")
write_results(results)
