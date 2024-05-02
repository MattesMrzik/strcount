#! /usr/bin/env python3
import sys
import pysam
import argparse
import pandas as pd


# TODO shouldn't the orientations be per line in config file?
def get_genome_str_graph(config=None,
                         reference_file=None,
                         repeat_orientation=None,
                         prefix_orientation=None,
                         suffix_orientation=None,
                         only_use_provided_fixes=False,
                         ucsc_browser_coords=False,
                         verbose=False):

    # read in the configs and store them in the respective variables
    configs = list()
    config_fh = open(config)
    header = config_fh.readline()

    for line in config_fh:
        if line.startswith("#"):
            continue
        if len(line) < 2:
            continue
        if verbose:
            sys.stderr.write("Parsing config file line: " + line.strip() + "\n")
        configs.append(line.rstrip().split()[0:7])

    config_fh.close()

    configs_dict = dict()

    # In this case I am extracting only 1 config but in case we have more than 1,
    # we can have a loop here to extract the other configs
    try:
        for chromosome, begin, end, name, repeat, prefix, suffix in configs:
            configs_dict[name] = [chromosome, int(begin), int(end), name, repeat, prefix, suffix]
    except Exception as e:
        sys.stderr.write("Error in config file: " + str(e) + "\n")
        sys.exit(1)

    # approach it chromosome wise and then have a variable for segments and one for
    # links for the chromosome and then in the end iterate over all and print.

    segments = dict()
    links = dict()

    c = 0

    # for column names see http://gfa-spec.github.io/GFA-spec/GFA1.html
    segments_cols = ["RecordType", "Name", "Sequence"]
    links_cols = ["RecordType", "From", "FromOrient", "To", "ToOrient", "Overlap"]
    for chr in pysam.FastxFile(reference_file):
        for key in configs_dict:
            config_chr, begin, end, name, repeat, prefix, suffix = configs_dict[key]
            if (chr.name == config_chr):
                # TODO what if the config file contains multiple entries for the same chromosome?
                if verbose:
                    sys.stderr.write("Processing: " + key + "\n")
                c = c + 1

                before_prefix_id = f"{chr.name}_before_prefix_{c}"
                prefix_id = f"{chr.name}_prefix_{c}"
                repeat_id = f"repeat_{c}"
                suffix_id = f"{chr.name}_suffix_{c}"
                after_suffix_id = f"{chr.name}_after_suffix_{c}"

                if not only_use_provided_fixes:

                    if ucsc_browser_coords:
                        begin = begin - 1
                        end = end
                    before_prefix_line = ["S", before_prefix_id, chr.sequence[: begin - len(prefix)]]
                    after_suffix_line = ["S", after_suffix_id, chr.sequence[end + len(suffix):]]
                prefix_line = ["S", prefix_id, prefix]
                repeat_line = ["S", repeat_id, repeat]
                suffix_line = ["S", suffix_id, suffix]

                dict_key = chr.name + "_" + str(c)

                if not only_use_provided_fixes:
                    segments[dict_key] = [before_prefix_line, prefix_line, repeat_line, suffix_line, after_suffix_line]
                else:
                    segments[dict_key] = [prefix_line, repeat_line, suffix_line]

                if verbose:
                    sys.stderr.write(f"Segments for chromosome {chr.name}\n")
                    segments_df = pd.DataFrame(segments[dict_key], columns=segments_cols)
                    sys.stderr.write(segments_df.to_string() + "\n")

                cigar = "0M"  # the "to" segment follows directly after the "from segment"

                if not only_use_provided_fixes:
                    before_prefix_line = ["L", before_prefix_id, "+", prefix_id, prefix_orientation, cigar]
                    after_suffix_line = ["L", suffix_id, suffix_orientation, after_suffix_id, "+", cigar]

                prefix_line = ["L", prefix_id, prefix_orientation, repeat_id, repeat_orientation, cigar]
                repeat_line = ["L", repeat_id, repeat_orientation, repeat_id, repeat_orientation, cigar]
                suffix_line = ["L", repeat_id, repeat_orientation, suffix_id, suffix_orientation, cigar]

                if not only_use_provided_fixes:
                    links[dict_key] = [before_prefix_line, prefix_line, repeat_line, suffix_line, after_suffix_line]
                else:
                    links[dict_key] = [prefix_line, repeat_line, suffix_line]

                if verbose:
                    sys.stderr.write(f"Links for chromosome {chr.name}\n")
                    segments_df = pd.DataFrame(links[dict_key], columns=links_cols)
                    sys.stderr.write(segments_df.to_string() + "\n")

            else:
                sys.stderr.write("Not processing: " + key + "\n")
                # TODO This should no tbe needed
                # segments[chr.name] = sep.join(["S", chr.name, chr.sequence])
    return segments, links


# TODO use to_csv to print the dataframe instead of stdout?
def print_genome_str_graph(segments: dict, links: dict):
    print("H\tVN:Z:a.0")
    for key in segments:
        for entry in segments[key]:
            print("\t".join(entry))

    for key in links:
        for entry in links[key]:
            print("\t".join(entry))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', help='the ref file', required=True)
    parser.add_argument('--config', help='the config file', required=True)
    parser.add_argument('--repeat_orientation', help='the orientation of the repeat string. + or -',
                        required=False, default="+")
    parser.add_argument('--prefix_orientation', help='the orientation of the prefix, + or -',
                        required=False, default="+")
    parser.add_argument('--suffix_orientation', help='the orientation of the suffix, + or -',
                        required=False, default="+")
    parser.add_argument('--verbose', action='store_true', help='verbose')
    parser.add_argument('--only_use_provided_fixes', action='store_true',
                        help='only use the provided suffixes and prefixes, not the whole reference sequence.')
    parser.add_argument('--ucsc_browser_coords', action='store_true',
                        help='use UCSC browser coordinates, i.e. start counting at 1 and end is inclusive.')

    args = parser.parse_args()

    reference_file = args.ref
    config = args.config
    repeat_orientation = args.repeat_orientation
    prefix_orientation = args.prefix_orientation
    suffix_orientation = args.suffix_orientation
    verbose = args.verbose
    only_use_provided_fixes = args.only_use_provided_fixes
    ucsc_browser_coords = args.ucsc_browser_coords

    print_genome_str_graph(*get_genome_str_graph(config, reference_file, repeat_orientation, prefix_orientation,
                                                 suffix_orientation, only_use_provided_fixes, ucsc_browser_coords, verbose))
