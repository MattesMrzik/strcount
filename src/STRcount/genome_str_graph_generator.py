#! /usr/bin/env python3
import sys
import pysam
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--ref', help='the ref file', required=True)
parser.add_argument('--config', help='the config file', required=True)
parser.add_argument('--repeat_orientation', help='the orientation of the repeat string. + or -', required=False, default="+")
parser.add_argument('--prefix_orientation', help='the orientation of the prefix, + or -', required=False, default="+")
parser.add_argument('--suffix_orientation', help='the orientation of the suffix, + or -', required=False, default="+")
parser.add_argument('--verbose', action='store_true', help='verbose')

args = parser.parse_args()

reference_file = args.ref
config = args.config
repeat_orientation = args.repeat_orientation
prefix_orientation = args.prefix_orientation
suffix_orientation = args.suffix_orientation
verbose = args.verbose

# read in the configs and store them in the respective variables
configs = list()
config_fh = open(config)
header = config_fh.readline()

for line in config_fh:
    configs.append(line.rstrip().split()[0:7])

configs_dict = dict()

# In this case I am extracting only 1 config but in case we have more than 1,
# we can have a loop here to extract the other configs
for chromosome, begin, end, name, repeat, prefix, suffix in configs:
    configs_dict[name] = [chromosome, begin, end, name, repeat, prefix, suffix]

print("H\tVN:Z:a.0")

# approach it chromosome wise and then have a variable for sides and one for
# links for the chromosome and then in the end iterate over all and print.

sides = dict()
links = dict()

c = 0

# for column names see http://gfa-spec.github.io/GFA-spec/GFA1.html
sides_cols = ["RecordType", "Name", "Sequence"]
links_cols = ["RecordType", "From", "FromOrient", "To", "ToOrient", "Overlap"]
sep = " ; "
for chr in pysam.FastxFile(reference_file):
    for key in configs_dict:
        if (chr.name == configs_dict[key][0]):
            # TODO what if the config file contains multiple entries for the same chromosome?
            if verbose:
                sys.stderr.write("Processing: " + key + "\n")
            c = c + 1

            before_prefix = ["S", chr.name + "_before_prefix", chr.sequence[: int(begin)-(len(prefix)+1)]]
            prefix = ["S", chr.name + "_prefix", chr.sequence[int(begin)-(len(prefix)+1): int(begin)-1]]
            repeat = ["S", "repeat_" + str(c), configs_dict[key][4]]
            suffix = ["S", chr.name + "_suffix", chr.sequence[int(end)+1: (int(end)+len(suffix)+1)]]
            after_suffix = ["S", chr.name + "_after_suffix", chr.sequence[int(end)+(len(suffix)+1):]]

            sides[chr.name] = [before_prefix, prefix, repeat, suffix, after_suffix]

            if verbose:
                sys.stderr.write(f"Sides for chromosome {chr.name}\n")
                sides_df = pd.DataFrame([before_prefix, prefix, repeat, suffix, after_suffix], columns=sides_cols)
                sys.stderr.write(sides_df.to_string() + "\n")

            before_prefix = ["L", chr.name + "_before_prefix", "+", chr.name + "_prefix", prefix_orientation, "*"]
            prefix = ["L", chr.name + "_prefix", prefix_orientation, "repeat_" + str(c), repeat_orientation, "*"]
            repeat = ["L", "repeat_" + str(c), repeat_orientation, "repeat_" + str(c), repeat_orientation, "*"]
            suffix = ["L", "repeat_" + str(c), repeat_orientation, chr.name + "_suffix", suffix_orientation, "*"]
            after_suffix = ["L", chr.name + "_suffix", suffix_orientation, chr.name + "_after_suffix", "+", "*"]

            links[chr.name] = [before_prefix, prefix, repeat, suffix, after_suffix]

            if verbose:
                sys.stderr.write(f"Links for chromosome {chr.name}\n")
                sides_df = pd.DataFrame([before_prefix, prefix, repeat, suffix, after_suffix], columns=links_cols)
                sys.stderr.write(sides_df.to_string() + "\n")

        else:
            sys.stderr.write("Not processing: " + key + "\n")
            # This should no tbe needed
            # sides[chr.name] = sep.join(["S", chr.name, chr.sequence])

# TODO use to_csv to print the dataframe instead of stdout?
for key in sides:
    for entry in sides[key]:
        print("\t".join(entry[1:]))


for key in links:
    for entry in links[key]:
        print("\t".join(entry[1:]))

if (c == 0):
    print("Error : The chromosome of the config file was not found in the reference, please double check your inputs.")
