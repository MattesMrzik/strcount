#! /usr/bin/env python

import argparse
import os
import logging

parser = argparse.ArgumentParser(description='Software tool to analyse STR loci from long read data. STRcount can count the number of repeats in a repeat expansion and give you the count in a tabular format for further downstream analysis.')
parser.add_argument('--reference', help='the reference from which the STR Graph will be generated', required=True)
parser.add_argument('--fastq', help='the baseaclled reads in fastq format', required=True)
parser.add_argument('--config', help='the config file', required=True)
parser.add_argument('--output', help='the output file', required=True)
parser.add_argument('--min-identity', type=float, default=0.50, help='only use reads with identity greater than this', required=False)
parser.add_argument('--min-aligned-fraction', type=float, default=0.8, help='require alignments cover this proportion of the query sequence', required=False)
parser.add_argument('--write-non-spanned', action='store_true', default=False, help='do not require the reads to span the prefix/suffix region', required=False)
parser.add_argument('--repeat_orientation', help='the orientation of the repeat string. + or -', required=False, default="+")
parser.add_argument('--prefix_orientation', help='the orientation of the prefix, + or -', required=False, default="+")
parser.add_argument('--suffix_orientation', help='the orientation of the suffix, + or -', required=False, default="+")
parser.add_argument('--cleanup', help='do you want to clean up the temporary file?', required=False, default="yes")
parser.add_argument('--output_directory', help='the output directory for all output and temporary files', required=False, default="./")
parser.add_argument('--multiseed-DP', help='Aligner option', required=False, default="")
parser.add_argument('--precise-clipping', help='Aligner option: use arg as the identity threshold for a valid alignment.', required=False, default="")
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for GraphAligner (default: 1)')
parser.add_argument('--verbose', action='store_true', help='verbose')
parser.add_argument('--only_use_provided_fixes', action='store_true', help='only use the provided suffixes and prefixes, not the complete reference sequence.')
parser.add_argument('--use_fixed_len_before_and_after_fixes', action='store_true', help='Use a fixed length for the before or after prefix or suffix sequences and not using the complete reference sequence.')

args = parser.parse_args()

ref = args.reference
reads = args.fastq
config_file = args.config
min_id = args.min_identity
min_align = args.min_aligned_fraction
span = args.write_non_spanned
rep_orientation = args.repeat_orientation
pre_orientation = args.prefix_orientation
suf_orientation = args.suffix_orientation
cleanup_flag = 0 if args.cleanup == "no" else 1
out_dir = args.output_directory
out_file = args.output
seed = args.multiseed_DP
id_thres = args.precise_clipping
threads = args.threads
verbose = args.verbose
only_use_provided_fixes = args.only_use_provided_fixes
use_fixed_len_before_and_after_fixes = args.use_fixed_len_before_and_after_fixes

if use_fixed_len_before_and_after_fixes:
    print("Error: --use_fixed_len_before_and_after_fixes is not implemented yet. Please remove this flag and try again.")
    exit(1)

#cwd = os.getcwd()

def print_red(s):
    print("\033[31m" + s + "\033[0m")

def main():
    if min_id:
        min_id_arg = f"--min-identity {min_id}"
    else:
        min_id_arg = ""

    if span:
        span_arg = f"--write-non-spanned {span}"
    else:
        span_arg = ""

    if min_align:
        min_align_arg = f"--min-aligned-fraction {min_align}"
    else:
        min_align_arg = ""

    if rep_orientation:
        rep_orientation_arg = f"--repeat_orientation {rep_orientation}"
    else:
        rep_orientation_arg = ""

    if pre_orientation:
        pre_orientation_arg = f"--prefix_orientation {pre_orientation}"
    else:
        pre_orientation_arg = ""

    if suf_orientation:
        suf_orientation_arg = f"--suffix_orientation {suf_orientation}"
    else:
        suf_orientation_arg = ""

    if seed:
        seed_arg = f"--multiseed-DP {seed}"
    else:
        seed_arg = ""

    if id_thres:
        id_thres_arg = f"--precise-clipping {id_thres}"
    else:
        id_thres_arg = ""

    if threads:
        threads_arg = f"-t {threads}"
    else:
        threads_arg = ""

    if verbose:
        verbose_arg = "--verbose"
    else:
        verbose_arg = ""

    if only_use_provided_fixes:
        only_use_provided_fixes_arg = "--only_use_provided_fixes"
    else:
        only_use_provided_fixes_arg = ""

    create_tmp_folder = os.system(f"mkdir {out_dir}tmp")

    optional_path_prefix = ""
    optional_path_prefix = "./src/STRcount/" # If you dont want to install to run the tool

    command = f"{optional_path_prefix}genome_str_graph_generator.py --ref {ref} --config {config_file} {rep_orientation_arg} {pre_orientation_arg} {suf_orientation_arg} {verbose_arg} {only_use_provided_fixes_arg} > {out_dir}tmp/genome_str_graph.gfa"
    print_red(command)
    str_graph_generator = os.system(command)

    if str_graph_generator == 0:
        logging.info("STR Reference Graph has been generated")
    else:
        logging.error("Error while Generating STR genome graph")

    command = f"GraphAligner {seed_arg} {id_thres_arg} -g {out_dir}tmp/genome_str_graph.gfa -f {reads} -a {out_dir}tmp/alignment.gaf -x vg {threads_arg} {verbose_arg}"
    print_red(command)
    ga_align = os.system(command)

    if ga_align == 0:
        logging.info("Reads aligned to Reference Graph")
    else:
        logging.error("Error in aligning the reads to the reference graph")

    command = f"{optional_path_prefix}parse_gaf.py --input {out_dir}tmp/alignment.gaf {min_id_arg} {min_align_arg} {span_arg} {verbose_arg} > {out_dir}{out_file}"
    print_red(command)
    parse_gaf = os.system(command)

    if parse_gaf == 0:
        logging.info("A read wise count has been generated")
    else:
        logging.error("Error in parsing the alignments")

if __name__ == "__main__":
    main()
