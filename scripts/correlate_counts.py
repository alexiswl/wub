#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import numpy as np
from scipy.stats import spearmanr
from scipy import stats
from collections import OrderedDict
import pandas as pd
from os import path
from wub.vis import report
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Correlate counts produced by multiple runs of bam_count_reads.py.""")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_multi_qc.pdf).", default="correlate_counts.pdf")
parser.add_argument(
    'counts', metavar='input_counts', nargs='*', type=str, help="Input counts as tab separated files.")


def load_counts(counts):
    """Load statistics from pickle files.

    :param pickles: List of count files.
    :returns: OrderedDict of count data frames per dataset.
    :rtype: OrderedDict
    """
    stats = OrderedDict()
    for count_file in counts:
        name = path.basename(count_file).rsplit('.', 1)[0]
        stats[name] = pd.read_csv(count_file, sep="\t")
    return stats


def _get_reference_set(dfs):
    """Get list of all references."""
    references = set()
    for df in six.itervalues(dfs):
        references = references.union(set(df['Reference']))
    return sorted(list(references))


def join_counts(counts):
    """Join count data frames.
    :param counts: Dicttionary of data frames.
    :returns: Mergeid data frame.
    :rtype: DataFrame
    """
    references = _get_reference_set(counts)
    res_dict = OrderedDict({'Reference': references})

    for dataset in six.iterkeys(counts):
        res_dict[dataset] = []

    for i, ref in enumerate(references):
        for dataset in six.iterkeys(counts):
            tmp = counts[dataset][counts[dataset]['Reference'] == ref]
            if len(tmp) == 0:
                tmp = 0
            elif len(tmp) == 1:
                tmp = int(tmp['Count'])
            else:
                raise Exception("Multiple rows for single reference in {}".format(dataset))
            res_dict[dataset].append(tmp)
    return pd.DataFrame(res_dict)


def _corrfunc(x, y, **kws):
    """ Annotate grid with correaltion coefficient.
    Solution from http://stackoverflow.com/a/30942817
    """
    r, _ = stats.spearmanr(x, y)
    ax = plotter.plt.gca()
    ax.annotate("R = {:.2f}".format(r),
                xy=(.1, .9), xycoords=ax.transAxes)


if __name__ == '__main__':
    args = parser.parse_args()
    plotter = report.Report(args.r)

    counts = load_counts(args.counts)
    joint_df = join_counts(counts)

    # Solution from http://stackoverflow.com/a/30942817
    g = sns.PairGrid(joint_df, palette=["red"])
    g.map_upper(plotter.plt.scatter, s=10)
    g.map_diag(sns.distplot, kde=False)
    g.map_lower(sns.kdeplot, cmap="Blues_d")
    g.map_lower(_corrfunc)
    g.map_upper(_corrfunc)

    plotter.plt.tight_layout()
    plotter.pages.savefig()

    plotter.close()
