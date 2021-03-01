[![Build Status](https://travis-ci.com/timodonnell/yabul.svg?branch=main)](https://travis-ci.com/timodonnell/yabul)
# yabul
Yet Another Bioinformatics Utilities Library

This is a small collection of Python functions for working with protein, DNA,
and RNA sequences. We use [pandas](https://pandas.pydata.org/) data frames
wherever possible. 

Yabul currently supports:
* Reading and writing FASTAs
* Pairwise local and global sequence alignment (uses [parasail](https://github.com/jeffdaily/parasail))

Requires Python 3.6+.
 
## Installation
Install using pip:

```
$ pip install yabul
```

You can run the unit from a checkout of the repo as follows:

```
$ pip install pytest
$ pytest
```

## Example

### Reading and writing FASTAs
The [read_fasta](https://yabul.readthedocs.io/en/latest/yabul/yabul.html#yabul.read_fasta)
function returns a [`pandas.DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html):
```
>>> import yabul
>>> df = yabul.read_fasta("test/data/cov2.fasta")
>>> df.head(3)
                                                             description                                           sequence
id
sp|P0DTC2|SPIKE_SARS2  sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein OS=Se...  MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSS...
sp|P0DTD1|R1AB_SARS2   sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab...  MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHL...
sp|P0DTC1|R1A_SARS2    sp|P0DTC1|R1A_SARS2 Replicase polyprotein 1a O...  MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHL...
```

The [write_fasta](https://yabul.readthedocs.io/en/latest/yabul/yabul.html#yabul.write_fasta) function takes 
(name, sequence) pairs:
```
>>> yabul.write_fasta("out.fasta", [("protein1", "TEST"), ("protein2", "HIHI")])
>>> yabul.write_fasta("out2.fasta", df.sequence.items())
```

### Sequence alignment
The [align_pair](https://yabul.readthedocs.io/en/latest/yabul/yabul.html#yabul.align_pair) function will give a local (Smith-Waterman) and global
(Needleman-Wunsch) alignment of two sequences. It returns a pandas.Series
with the aligned sequences.

By default, the alignment is global:
```
>>> yabul.align_pair("AATESTDD", "TEST")
query             AATESTDD
reference         --TEST--
correspondence      ||||
score                   -5
dtype: object
```

To do a local alignment, pass `local=True`.
```
>>> yabul.align_pair("AATESTDD", "TEST", local=True)
query             TEST
reference         TEST
correspondence    ||||
score               19
dtype: object
```

## Dependencies
The alignment routine is a thin wrapper around the Smith-Waterman and
Needleman-Wunsch implementations from [parasail](https://github.com/jeffdaily/parasail).

## Contributing
We welcome contributions of well-documented code to read and write common
bioinformatics file formats using pandas objects. Please include unit tests
in your PR. Additional functionality like multiple sequence alignment would
also be nice to add.

## Releasing
To push a new release to PyPI:
* Make sure the package version specified in [`__init__.py`](https://github.com/timodonnell/yabul/blob/main/yabul/__init__.py)
is a new version greater than what's on [PyPI](https://pypi.org/project/yabul/).
* Tag a new release on GitHub matching this version

Travis should deploy the release to PyPI automatically.

Documentation at https://yabul.readthedocs.io/en/latest/ should update automatically on commit.

To build the documentation locally, run:

```
$ cd docs
$ pip install -r requirements.txt
$ sphinx-build -b html . _build
```

