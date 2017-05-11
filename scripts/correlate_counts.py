#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import sys
import numpy as np
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
    '-c', metavar='corr_type', type=str, help="Correlation statistic - spearman or pearson (spearman).", default="spearman")
parser.add_argument(
    '-L', action="store_true", help="Log transform data.", default=False)
parser.add_argument(
    'counts', metavar='input_counts', nargs='*', type=str, help="Input counts as tab separated files.")


def load_counts(counts, log_transform):
    """Load statistics from pickle files.

    :param pickles: List of count files.
    :returns: OrderedDict of count data frames per dataset.
    :rtype: OrderedDict
    """
    stats = OrderedDict()
    for count_file in counts:
        name = path.basename(count_file).rsplit('.', 1)[0]
        if log_transform:
            name = 'log(' + name + '+1)'
        stats[name] = pd.read_csv(count_file, sep="\t")
        stats[name]['Count'] = np.log(stats[name]['Count'] + 1)
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
    if args.c == 'spearman':
        r, _ = stats.spearmanr(x, y)
        corr_type = 'Rho'
    elif args.c == 'pearson':
        r, _ = stats.pearsonr(x, y)
        corr_type = 'r'
    else:
        raise Exception('Invalid correlation statistic.')
    correlations.append(r)
    ax = plotter.plt.gca()
    ax.annotate("{} = {:.2f}".format(corr_type, r),
                xy=(.1, .9), xycoords=ax.transAxes)


if __name__ == '__main__':
    args = parser.parse_args()
    plotter = report.Report(args.r)

    if len(args.counts) == 0:
        sys.stderr.write("No count files given!\n")
        sys.exit(1)

    counts = load_counts(args.counts, args.L)
    joint_df = join_counts(counts)
    correlations = []

    # Solution from http://stackoverflow.com/a/30942817
    g = sns.PairGrid(joint_df, palette=["red"])
    g.map_upper(plotter.plt.scatter, s=10)
    g.map_diag(sns.distplot, kde=False)
    g.map_lower(sns.kdeplot, cmap="Blues_d")
    g.map_lower(_corrfunc)
    g.map_upper(_corrfunc)
    plotter.plt.tight_layout()
    plotter.pages.savefig()

    plotter.plt.clf()
    correlations = pd.DataFrame(
        {"Distribution of correlation coefficients": correlations})
    sns.boxplot(correlations)
    plotter.plt.tight_layout()
    plotter.pages.savefig()

    plotter.close()
