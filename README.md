# Virconsens
Tool to create a consensus sequence from mapped virus Nanopore data

## Installation

1. Clone this repository and ``cd Virconsens``
2. ``conda env create -f environment.yml --name virconsens``
3. ``conda activate virconsens``
4. ``pip install .``

## Full usage

```
usage: virconsens.py [-h] -b BAM -o OUT -n OUTNAME -r REFERENCE
                     [-vf VARIANTFILE] [-c CORES] [-d MINDEPTH] [-af MINAF]

Virconsens

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     BAM file from which to create a consensus
  -o OUT, --out OUT     Output path for consensus fasta
  -n OUTNAME, --outname OUTNAME
                        Name to be given to the output consensus sequence
  -r REFERENCE, --reference REFERENCE
                        Reference genome fasta file
  -vf VARIANTFILE, --variantfile VARIANTFILE
                        Output path for variant tsv file
  -c CORES, --cores CORES
                        Number of cores to use for processing
  -d MINDEPTH, --mindepth MINDEPTH
                        Minimal depth at which to not consider any alternative
                        alleles
  -af MINAF, --minAF MINAF
                        Minimal allele frequency to output
```