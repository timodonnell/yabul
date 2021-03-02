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
    # We (mis-) use pandas to parse the file.
    df = pandas.read_csv(
        filename,
        header=None,
        skip_blank_lines=True,
        quoting=3,  # QUOTE_NONE
        comment=';',  # Fasta comment lines start with ';'
        sep="\0",  # null separator: never split, always read one column
    )
    assert df.shape[1] == 1  # one column
    df.columns = ["raw"]

    is_header_line = df.raw.str.startswith(">")

    # Assign a new index for every line starting with ">"
    df.loc[is_header_line, "idx"] = 1
    df.idx = df.idx.cumsum()
    df.idx = df.idx.fillna(method="ffill").astype(int)

    idx_to_description = df.loc[is_header_line].set_index("idx").raw.str.slice(1)

    df.loc[~is_header_line, "sequence_piece"] = df.raw
    df.sequence_piece = df.sequence_piece.fillna("")

    result = df.groupby("idx").sequence_piece.apply("".join).to_frame()
    result.columns = ["sequence"]
    result.insert(0, "description", idx_to_description)
    result.index = result.description.str.split().str.get(0)
    result.index.name = "id"
    return result
