# -*- coding: utf-8 -*-
"""Yet uncategorised utility functions."""

import pickle
import os.path
import sys
from itertools import chain


def get_fname(fname):
    """ get the file name without extension.

    :param fname: file name
    :return: file name
    :rtype: str

    """
    return os.path.splitext(os.path.basename(fname))[0]


def get_extension(fname):
    """ get the file extension.

    :param fname: file name
    :return: file extention
    :rtype: str format '.*'

    """
    return os.path.splitext(os.path.basename(fname))[1]


def _getextension(fast):
    """ finds and check for the correct extension. If extension is not correct it will return Exception and exit.

    :param fast: fastq or fasta file
    :return: "fastq" or "fasta"
    :rtype: str

    """

    extension = get_extension(fast)
    if extension in ('.fa', '.fasta'):
        extension = "fasta"
    elif extension in ('.fq', '.fastq'):
        extension = "fastq"
    else:
        raise Exception('Incorrect file format')
        exit()
        # print >> sys.stderr, "Incorrect file format"
    return extension


def mkdir(path):
    """ if the dir does not exists it create it

    :param path: dir path
    :return: path
    :rtype: str

    """
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def pickle_load(fname):
    """ Load object from pickle.

    :param fname: Input pickle file name.
    :returns: Object loaded from pickle file.
    :rtype: object

    """
    fh = open(fname, 'rb')
    data = pickle.load(fh)
    fh.close()
    return data


def pickle_dump(obj, fname):
    """Pickle object to file.

    :param obj: Object to be pickled.
    :fname: Output file name.
    :returns: The name of output file.
    :rtype: str

    """
    fh = open(fname, 'wb')
    pickle.dump(obj, fh)
    fh.flush()
    fh.close()
    return fname


def merge_pickles(obj_list):
    """

    :param obj_list: list of imported pickle objects
    :return: merged_pickles
    """

    # Ensure all pickles have the same tag
    if len(set([pickle['tag'] for pickle in obj_list])) > 1:
        sys.exit("Multiple tags not yet supported")

    # Create new object
    merged = dict()

    # Set tag
    merged['tag'] = obj_list[0]['tag']

    # Set read stats
    # Standard set of keys for each object
    read_stats_keys = list(obj_list[0]['read_stats'].keys())
    # Initialise object
    merged['read_stats'] = {}
    # Iterate through each key for read_stats
    for key in read_stats_keys:
        # The following attributes are summed
        if key in ['mapped', 'unmapped']:
            merged['read_stats'][key] = sum([obj['read_stats'][key] for obj in obj_list])
        # The following attributes are appended
        else:
            merged['read_stats'][key] = list(chain(*(obj['read_stats'][key] for obj in obj_list)))

    # Set error stats
    merged['error_stats'] = {}
    # First get all keys (from AAA to NNN)
    all_error_keys = list(set(chain(*(list(obj['error_stats'].keys()) for obj in obj_list))))
    # Each AAA is a dict (A, C, G, T, '-' *)
    all_nucs = ['A', 'C', 'G', 'T', '-', '*']
    for error_key in all_error_keys:
        merged['error_stats'][error_key] = {}
        for nuc in all_nucs:
            # Get values
            nuc_value = sum(obj['error_stats'][error_key][nuc]
                            for obj in obj_list
                            if error_key in obj['error_stats'].keys()
                            and nuc in obj['error_stats'][error_key].keys())
            if not nuc_value == 0:
                merged['error_stats'][error_key][nuc] = nuc_value

    # Set indel_stats:
    # Standard for each stats
    indel_stats_keys = ['insertion_lengths', 'deletion_lengths', 'insertion_composition']
    all_nucs = ['A', 'C', 'G', 'T']
    merged['indel_stats'] = {}
    for indel_key in indel_stats_keys:
        merged['indel_stats'][indel_key] = {}
        if indel_key == 'insertion_composition':
            for nuc in all_nucs:
                merged['indel_stats'][indel_key][nuc] = sum(obj_list['indel_stats'][indel_key][nuc]
                                                            for obj in obj_list)
        else:
            # Get all values for this key
            all_insertion_lengths = list(set(chain(*(obj['indel_stats'][indel_key].keys() for obj in obj_list))))
            # Sum for each insertion length
            for length in all_insertion_lengths:
                merged['indel_stats'][indel_key][length] = sum(obj['indel_stats'][indel_key][length]
                                                               for obj in obj_list
                                                               if length in obj['indel_stats'][indel_key].keys())

    # Set the base stats
    # Standard set of keys
    base_stats_keys = list(obj_list[0]['base_stats'].keys())
    merged['base_stats'] = {}
    for key in base_stats_keys:
        if key not in ['identity', 'accuracy']:
            merged['base_stats'][key] = sum(obj['base_stats'][key] for obj in obj_list)
    # Now calculate identity and accuracy stats
    merged['base_stats']['identity'] = sum(merged['base_stats']['match']) / \
                                       (sum(merged['base_stats']['match']) +
                                        sum(merged['base_stats']['mismatch']))
    merged['base_stats']['accuracy'] = sum(merged['base_stats']['match']) / \
                                       (sum(merged['base_stats']['match']) +
                                        sum(merged['base_stats']['mismatch']) +
                                        sum(merged['base_stats']['insertion']) +
                                        sum(merged['base_stats']['deletion']))

    return merged
