#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from wub.mappers import lastal

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare a set of reference sequences (genome) to another set (target assembly) using lastal alignment.
    Accuracy is the total number of matched bases divided by total alignment length. Coverage is total reference covered
    by alignment divided by total length of reference.

    Caveats:
     - The lastal alignments are filtered so only the best scoring alignment is kept per query. Hence some shorter valid
     alignments might be discarded causing an underestimation of coverage.
    """)
parser.add_argument('ref', metavar='reference_fasta', type=str, help="Reference fasta.")
parser.add_argument('target', metavar='target_fasta', type=str, help="Target fasta.")


if __name__ == '__main__':
    args = parser.parse_args()

    stats = lastal.compare_genomes_lastal(args.ref, args.target, {}, cleanup=True)

    global_accuracy = (stats['aln_length'].sum() - stats['substitutions'].sum() -
                       stats['deletions'].sum() - stats['insertions'].sum()) / float(stats['aln_length'].sum())
    global_coverage = stats['ref_aln_len'].sum() / float(stats['ref_len'].sum())

    sys.stdout.write("Accuracy\tCoverage\n")
    sys.stdout.write("{}\t{}\n".format(global_accuracy, global_coverage))