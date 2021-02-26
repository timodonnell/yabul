# yabul
Yet Another Bioinformatics Utilities Library

This is a small collection of Python functions for working with protein, DNA,
and RNA sequences. We use [pandas](https://pandas.pydata.org/) data frames
wherever possible.

Yabul currently supports:
* Pairwise local and global sequence alignment (uses [parasail](https://github.com/jeffdaily/parasail))
* Reading and writing FASTAs
 
## Installation
From a checkout run:

```
$ pip install .
```

You can run the unit tests with:

```
$ pip install pytest
$ pytest test/
```

## Example

The [align_pair] function will give a local (Smith-Waterman) and global
(Needleman-Wunsch) alignment of two sequences. It returns a pandas.Series
with the aligned sequences.
```
>>> yabul.align_pair("AATESTDD", "TEST")
query             AATESTDD
reference         --TEST--
correspondence      ||||
score                   -5
dtype: object
```

Local alignment:
```
>>> yabul.align_pair("AATESTDD", "TEST", local=True)
query             TEST
reference         TEST
correspondence    ||||
score               19
dtype: object
```


The [read_fasta] function returns a pandas.DataFrame:
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

The [write_fasta]() function takes (name, sequence) pairs:
```
>>> yabul.write_fasta("out.fasta", [("protein1", "TEST"), ("protein2", "HIHI")])
```

## Dependencies
The alignment routine is a thin wrapper around the Smith-Waterman and
Needleman-Wunsch implementations from [parasail](https://github.com/jeffdaily/parasail).
