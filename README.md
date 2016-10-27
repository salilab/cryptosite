[![Build Status](https://travis-ci.org/salilab/cryptosite.svg?branch=master)](https://travis-ci.org/salilab/cryptosite)
[![codecov](https://codecov.io/gh/salilab/cryptosite/branch/master/graph/badge.svg)](https://codecov.io/gh/salilab/cryptosite)

CryptoSite is a computational tool for predicting the location of cryptic
binding sites in proteins and protein complexes.

A [web interface](https://salilab.org/cryptosite/) is also available.

# Usage

Input files:

- param.txt:
  job name
  email address
  chain

- input.pdb:
  input PDB file

- input.seq:
  input sequence in fasta format or just sequence


Run by executing this line:

  /salilab/diva1/home/modeller/modpy.sh python mainer.py
