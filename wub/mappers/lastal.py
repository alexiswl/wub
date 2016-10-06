import os
import re
from collections import namedtuple
import subprocess
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO

lastdb_suffixes = ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']

record_template = 'score r_name r_start r_aln_len r_strand r_len r_aln '
record_template += 'q_name q_start q_aln_len q_strand q_len q_aln'
LastRecord = namedtuple('LastRecord', record_template)


def lastdb(ref_dir, ref_name, ref, executable='lastdb', **kwargs):
    """
    Runs lastdb on ref within ref_dir using the label ref_name if
    any errors thrown during runtime, files are checked for existence
    if all files accounted for, successful=False but no errors thrown.
    Otherwise, IOError or CalledProcessError thrown.


    :param ref_dir: directory you will find lastdb files in
    :param ref_name: name of the lastdb files e.g. a for a.prj..
    :param ref: filepath for reference file
    :param executable: path/executable for lastdb e.g. ont_lastdb
    :param kwargs: any -[arg] wanted see lastdb -h for details
    :return: True/False is successful with no errors and command run
    :raises: `IOError` if files don't exist
    :raises: `subprocess.CalledProcessError` for errors during runtime
    """
    # taken from metrichor-bio/alignment:

    ref_dir = os.path.abspath(ref_dir)
    ref = os.path.abspath(ref)
    chdir = 'cd {}'.format(ref_dir)
    kwargs = ['-{} {}'.format(a, b) for a, b in kwargs.items()]
    kwargs_str = ' '.join(kwargs)
    cmd = '{} {} {} {}'.format(executable, ref_name, ref, kwargs_str)
    command = '{}; {}'.format(chdir, cmd)

    # Because a cd is done in subprocess, if that fails cwd is used.
    if not os.path.exists(ref_dir):
        raise IOError('Directory not found: {}'.format(ref_dir))

    try:
        _ = subprocess.check_output(command, shell=True,
                                    stderr=subprocess.STDOUT)
        successful = True
    except subprocess.CalledProcessError as e:
        successful = False
        if not os.path.exists(ref):
            raise IOError('Reference not found: {}'.format(ref))
        elif not utilities.find_executable(executable):
            raise IOError('Executable not found: {}'.format(executable))
        elif self.check_lastdb_file(ref_dir, ref_name):
            raise e

    return os.path.join(ref_dir, ref_name)


def check_lastdb_files(ref_dir, name):
    """
    Check that all lastdb files with `name` label exist within directory

    :param ref_dir: directory to check for lastdb files
    :param name: label to search for e.g. 'a' for a.prj
    :return: list of missing extensions, [] if none missing
    """
    # taken from metrichor-bio/alignment:
    files = 0
    missing = []
    for suffix in lastdb_suffixes:
        lastdb_file = '.'.join([name, suffix])
        lastdb_path = os.path.join(ref_dir, lastdb_file)
        files += os.path.exists(lastdb_path)
        if not os.path.exists(lastdb_path):
            missing.append(suffix)
    return missing


def lastal_align(database, query, executable='lastal', **kwargs):
    """
    Runs lastal via subprocess

    :param database: database prefix
    :param query: filepath for the query file
    :param kwargs: -[args] wanted for lastal e.g. v='' for verbosity
    :return: alignment output
    """
    # some code taken from metrichor-bio/alignment:
    ref_dir, ref_name = os.path.split(database)
    chdir_cmd = 'cd {}'.format(ref_dir)
    kwargs = ['-{} {}'.format(a, b) for a, b in kwargs.items()]
    kwargs_str = ' '.join(kwargs)

    cmd = ' '.join([executable, kwargs_str, ref_name, query])
    command = '{}; {}'.format(chdir_cmd, cmd)

    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

    # FIXME: reading whole alignment in memory, might not always be a good idea:
    return p.stdout.read()


def parse_lastal(res):
    """
    Parse raw lastal output records.
    :param res: Raw lastal results.
    :returns: Generator of lastal alignment records.
    :rtype: generator
    """
    lines = [l for l in res.split('\n') if not l.startswith('#')]
    re_score = re.compile('\s|=')
    re_aln = re.compile('\s+')
    for i, line in enumerate(lines):
        if line.startswith('a '):
            # Parse out score:
            score = int(re_score.split(line)[2])
            # Parse reference record:
            tmp = re_aln.split(lines[i + 1])[1:7]
            ref_name = tmp[0]
            ref_start = int(tmp[1])
            ref_aln_len = int(tmp[2])
            ref_strand = tmp[3]
            ref_len = int(tmp[4])
            ref_aln = tmp[5]
            # Parse query record:
            tmp = re_aln.split(lines[i + 2])[1:7]
            q_name = tmp[0]
            q_start = int(tmp[1])
            q_aln_len = int(tmp[2])
            q_strand = tmp[3]
            q_len = int(tmp[4])
            q_aln = tmp[5]
            # Yield record:
            record = LastRecord(score, ref_name, ref_start, ref_aln_len, ref_strand, ref_len,
                                ref_aln, q_name, q_start, q_aln_len, q_strand, q_len, q_aln)
            yield record