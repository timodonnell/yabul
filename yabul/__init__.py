"""
Yet Another Bioinformatics Utility Library
"""
__version__ = "0.0.2"
from .fasta import read_fasta, write_fasta
from .align import align_pair

__all__ = [
    "read_fasta",
    "write_fasta",
    "align_pair",
]
