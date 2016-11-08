[![docs](https://readthedocs.org/projects/cryptosite/badge/)](https://cryptosite.readthedocs.org/)
[![Build Status](https://travis-ci.org/salilab/cryptosite.svg?branch=master)](https://travis-ci.org/salilab/cryptosite)
[![codecov](https://codecov.io/gh/salilab/cryptosite/branch/master/graph/badge.svg)](https://codecov.io/gh/salilab/cryptosite)

CryptoSite is a computational tool for predicting the location of cryptic
binding sites in proteins and protein complexes.

A [web interface](https://salilab.org/cryptosite/) is also available.

# Dependencies

The CryptoSite tools expect to be able to find the following tools in standard
system paths (e.g. `PATH`, `PYTHONPATH`):

- [MODELLER](https://salilab.org/modeller/) plus the
  [SOAP-Protein library](https://salilab.org/SOAP/)
- [MUSCLE](http://www.drive5.com/muscle/)
- [DSSP](http://swift.cmbi.ru.nl/gv/dssp/)
- [fpocket (version 2)](http://fpocket.sourceforge.net/)
- [PatchDock](http://bioinfo3d.cs.tau.ac.il/PatchDock/)
- [IMP](https://integrativemodeling.org/)
- [NCBI BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  plus a local copy of the [UniProt database](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/)
  for it to search against
  (the local database should be a concatenation of the `uniprot_trembl.fasta`
  and `uniprot_sprot.fasta` files available from the UniProt website)
- [USEARCH](http://drive5.com/usearch/download.html)
- [Biopython](http://biopython.org/)
- [scikit-learn](http://scikit-learn.org/)

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
