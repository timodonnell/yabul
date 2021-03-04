"""
FASTA reading and writing
"""

from __future__ import print_function, division, absolute_import

import gzip
import textwrap

import pandas


def write_fasta(filename, sequences):
    """
    Write sequences to a FASTA.

    Parameters
    ----------
    filename : string
        File to write. If it ends with '.gz' the file will be gzip compressed.

    sequences : iterable of (name, sequence) pairs
        Sequences to write. Both name and sequence should be strings.

    """
    # wt is "open for writing in text mode"
    handle = (
        gzip.open(filename, "wt")
        if filename.endswith(".gz")
        else open(filename, "wt"))
    try:
        for name, sequence in sequences:
            handle.write(">")
            handle.write(name)
            handle.write("\n")
            for line in textwrap.wrap(sequence):
                handle.write(line)
                handle.write("\n")
            handle.write("\n")
    finally:
        handle.close()


def read_fasta(filename):
    """
    Parse a fasta file to a pandas DataFrame.

    Compression is supported (via pandas read_csv) and is inferred by
    extension: '.gz', '.bz2', '.zip', or '.xz'.

    Parameters
    ----------
    filename : string

    Returns
    -------
    pandas.DataFrame with columns "description" and "sequence". The index of the
    DataFrame is the "sequence ID", i.e. the first space-separated token of the
    description.
    """
    # We (mis-) use pandas to read the file.
    lines = pandas.read_csv(
        filename,
        header=None,
        skip_blank_lines=True,
        dtype=str,
        na_filter=False,
        quoting=3,  # QUOTE_NONE
        comment=';',  # Fasta comment lines start with ';'
        sep="\0",  # null separator: never split, always read one column
    )
    assert lines.shape[1] == 1  # one column
    lines.columns = ["raw"]

    result = []
    current_id = None
    pieces = []
    for line in lines.raw:
        if line.startswith(">"):
            if current_id is not None:
                result.append((current_id, "".join(pieces)))
                pieces.clear()
            current_id = line[1:]
        else:
            pieces.append(line)

    # Handle last sequence
    if current_id is not None:
        result.append((current_id, "".join(pieces)))

    result = pandas.DataFrame(result, columns=["description", "sequence"])
    result.index = result.description.str.split().str.get(0)
    result.index.name = "id"

    return result
