"""
FASTA reading and writing
"""

from __future__ import print_function, division, absolute_import

import gzip
import bz2
import textwrap
import os

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


def read_fasta(
        filename,
        max_sequences=None,
        chunksize=10000,
        progress=True,
        encoding='ascii'):
    """
    Parse a fasta file to a pandas DataFrame.

    Compression is supported and is inferred by extension: '.gz' or '.bz2'

    Parameters
    ----------
    filename : string
        Path to file to read
    chunksize : int, optional
        Number of lines to read at once
    progress : bool, optional
        If True and the tqdm package is available, a progress bar will be shown.
        The progrses bar shows progress through the file (in bytes) as it is read.
    encoding : string, optional
        Character encoding

    Returns
    -------
    pandas.DataFrame with columns "description" and "sequence". The index of the
    DataFrame is the "sequence ID", i.e. the first space-separated token of the
    description.
    """
  
    pbar = None
    raw_fd = None
    try:
        raw_fd = open(filename, "rb")
        compression = {
            "gz": "gzip",
            "bz2": "bz2",
        }.get(filename.split(".")[-1])

        if compression == "gzip":
            fd = gzip.open(raw_fd)
        elif compression == "bz2":
            fd = bz2.open(raw_fd)
        else:
            fd = raw_fd

        try:
            import tqdm
        except ImportError:
            progress = False

        if progress:
            # Monitor progress through file for progress indicator
            size = os.path.getsize(filename)
            pbar = tqdm.tqdm(total=size, unit='B', unit_scale=True, unit_divisor=1024)
            current_max = [0]  # list to make it mutable
            def progress_update(total):
                tell = raw_fd.tell()
                if tell > current_max[0]:
                    pbar.update(tell - current_max[0])
                    current_max[0] = tell
                    pbar.set_postfix({'seqs': total})
        else:
            def progress_update(total):
                pass

        # We (mis-) use pandas to read the file.
        iterator = pandas.read_csv(
            fd,
            header=None,
            skip_blank_lines=True,
            dtype=str,
            na_filter=False,
            chunksize=chunksize,
            encoding=encoding,
            encoding_errors='replace',
            quoting=3,  # QUOTE_NONE
            comment=';',  # Fasta comment lines start with ';'
            sep="\0",  # null separator: never split, always read one column
        )

        result = []
        current_id = None
        pieces = []
        dfs = []
        total = 0

        for lines in iterator:
            assert lines.shape[1] == 1  # one column
            lines.columns = ["raw"]

            for line in lines.raw:
                if line.startswith(">"):
                    if current_id is not None:
                        result.append((current_id, "".join(pieces)))
                        pieces.clear()
                        if max_sequences is not None and len(result) + total >= max_sequences:
                            break
                    current_id = line[1:]
                else:
                    pieces.append(line)
            
            dfs.append(pandas.DataFrame(result, columns=["description", "sequence"]))
            total += len(result)
            result.clear()
            progress_update(total)
            
            if max_sequences is not None and total >= max_sequences:
                break

    finally:
        if raw_fd is not None:
            raw_fd.close()
        if pbar is not None:
            pbar.close()

    # Handle last sequence
    if current_id is not None and total != max_sequences:
        result.append((current_id, "".join(pieces)))
        dfs.append(pandas.DataFrame(result, columns=["description", "sequence"]))
        total += len(result)
        result.clear()

    result = pandas.concat(dfs, ignore_index=True)
    result.index = result.description.str.split().str.get(0)
    result.index.name = "id"

    return result
