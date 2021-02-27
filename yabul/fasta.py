"""
Adapted from mhcflurry via pyensembl, github.com/openvax/pyensembl
Original pyensembl implementation by Alex Rubinsteyn.
"""

from __future__ import print_function, division, absolute_import

import gzip
import logging
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

    Parameters
    ----------
    filename : string
        If the filename ends with '.gz' is will be read as gzip'd input.

    Returns
    -------
    pandas.DataFrame with columns "description" and "sequence". The index of the
    DataFrame is the "sequence ID", i.e. the first space-separated token of the
    description.
    """
    reader = FastaParser()
    rows = reader.iterate_over_file(filename)
    result = pandas.DataFrame(
        rows,
        columns=["description", "sequence"])
    result.index = result.description.str.split().str.get(0)
    result.index.name = "id"
    return result


class FastaParser(object):
    """
    FastaParser object consumes lines of a FASTA file incrementally.
    """
    def __init__(self):
        self.current_id = None
        self.current_lines = []

    def iterate_over_file(self, fasta_path):
        """
        Generator that yields identifiers paired with sequences.
        """
        with self.open_file(fasta_path) as f:
            for line in f:
                line = line.rstrip()

                if len(line) == 0:
                    continue

                # have to slice into a bytes object or else get a single integer
                first_char = line[0:1]

                if first_char == b">":
                    previous_entry = self._current_entry()
                    self.current_id = self._parse_header_id(line)

                    if len(self.current_id) == 0:
                        logging.warning(
                            "Unable to parse ID from header line: %s", line)

                    self.current_lines = []

                    if previous_entry is not None:
                        yield previous_entry

                elif first_char == b";":
                    # semicolon are comment characters
                    continue
                else:
                    self.current_lines.append(line)

        # the last sequence is still in the lines buffer after we're done with
        # the file so make sure to yield it
        id_and_seq = self._current_entry()
        if id_and_seq is not None:
            yield id_and_seq

    def _current_entry(self):
        # when we hit a new entry, if this isn't the first
        # entry of the file then put the last one in the dictionary
        if self.current_id:
            if len(self.current_lines) == 0:
                logging.warning("No sequence data for '%s'", self.current_id)
            else:
                sequence = b"".join(self.current_lines).decode("ascii")
                return self.current_id, sequence

    @staticmethod
    def open_file(fasta_path):
        """
        Open either a text file or compressed gzip file as a stream of bytes.
        """
        if fasta_path.endswith("gz") or fasta_path.endswith("gzip"):
            return gzip.open(fasta_path, 'rb')
        else:
            return open(fasta_path, 'rb')

    @staticmethod
    def _parse_header_id(line):
        """
        Pull the transcript or protein identifier from the header line
        which starts with '>'
        """
        if len(line) <= 1:
            raise ValueError("No identifier on FASTA line")

        return line[1:].decode("ascii")