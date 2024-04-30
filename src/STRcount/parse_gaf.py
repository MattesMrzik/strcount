#! /usr/bin/env python

import argparse
import sys


class GraphAlignment:
    def __init__(self, read_name, strand, spanned, count, alignment_score, identity, aligned_fraction):
        self.read_name = read_name
        self.strand = strand
        self.spanned = spanned
        self.count = count
        self.alignment_score = alignment_score
        self.identity = identity
        self.aligned_fraction = aligned_fraction

    def print(self):
        sys.stderr.write(
            f"read_name: {self.read_name}, strand: {self.strand}, spanned: {self.spanned}, count: {self.count}, alignment_score: {self.alignment_score}, identity: {self.identity}, aligned_fraction: {self.aligned_fraction}\n")


parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from graphaligner', required=True)
parser.add_argument('--min-identity', type=float, default=0.50,
                    help='only use reads with identity greater than this', required=False)
parser.add_argument('--min-aligned-fraction', type=float, default=0.8,
                    help='require alignments cover this proportion of the query sequence', required=False)
parser.add_argument('--write-non-spanned', action='store_true', default=False,
                    help='do not require the reads to span the prefix/suffix region', required=False)
parser.add_argument('--verbose', action='store_true', help='verbose')

args = parser.parse_args()

input_file = open(args.input)

alignments = dict()
with open(args.input) as f:
    for record in f:
        fields = record.rstrip().split("\t")

        if args.verbose:
            sys.stderr.write(f"Processing alignment: {', '.join(fields)}\n")

        read_id = fields[0].split(' ')[0]  # remove FASTQ metadata that GraphAligner emits

        tags = dict()
        for t in fields[12:]:
            (key, data_type, value) = t.split(":")
            tags[key] = value

        align_score = float(tags["AS"])
        identity = float(tags["id"])

        query_len = int(fields[1])
        query_start = int(fields[2])
        query_end = int(fields[3])
        query_af = float(query_end - query_start) / float(query_len)

        path_str = fields[5]
        path_dir = path_str[0]
        segments = path_str.split(path_dir)
        has_prefix = False
        has_suffix = False

        count = 0
        for s in segments:
            if s.find("prefix") >= 0:
                has_prefix = True
            if s.find("suffix") >= 0:
                has_suffix = True
            count += s.find("repeat") >= 0

        valid = has_prefix and has_suffix
        strand = fields[4]
        if args.verbose:
            if strand != "+":
                sys.stderr.write("Strand is not \"+\", exiting.\n")
        assert (strand == "+")
        if path_dir == "<":
            strand = "-"
        ga = GraphAlignment(read_id, strand, valid, count, align_score, identity, query_af)
        if args.verbose:
            sys.stderr.write("Found alignment: \n")
            ga.print()
        if not valid and not args.write_non_spanned:
            if args.verbose:
                sys.stderr.write("Skipping alignment because it does not span the prefix/suffix region\n")
            continue
        if ga.identity < args.min_identity:
            if args.verbose:
                sys.stderr.write("Skipping alignment because identity is too low\n")
            continue
        if ga.aligned_fraction < args.min_aligned_fraction:
            if args.verbose:
                sys.stderr.write("Skipping alignment because aligned fraction is too low\n")
            continue
        if ga.read_name not in alignments or alignments[ga.read_name].alignment_score < ga.alignment_score:
            if args.verbose:
                sys.stderr.write("Keeping alignment bc it has better score than previous alignment\n")
            alignments[ga.read_name] = ga

print("\t".join(["read_name", "strand", "spanned", "count", "align_score", "identity", "query_aligned_fraction"]))
for ga in alignments.values():
    print("%s\t%s\t%d\t%d\t%.1f\t%.3f\t%.3f" % (ga.read_name, ga.strand, ga.spanned,
          ga.count, ga.alignment_score, ga.identity, ga.aligned_fraction))
