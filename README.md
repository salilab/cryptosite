Dependencies:
=============

1. blast	- http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
2. muscle	- http://www.drive5.com/muscle/
3. usearch	- http://drive5.com/usearch/
4. CHASA.py	- http://folding.chemistry.msstate.edu/chasa/chasa_file.py
5. Fpocket	- http://fpocket.sourceforge.net/
6. NumPy	- http://www.numpy.org/
7. SciPy	- http://www.scipy.org/
8. SciKit-Learn	- http://scikit-learn.org/stable/
9. BioPython	- http://biopython.org/wiki/Main_Page
10. Modeller	- https://salilab.org/modeller/
11. DSSP	- http://swift.cmbi.ru.nl/gv/dssp/
12. Uniprot protein sequence database	- http://www.uniprot.org/downloads
13. AllosMod	- https://github.com/salilab/allosmod
14. PatchDock	- http://bioinfo3d.cs.tau.ac.il/PatchDock/

Input files:
============

- param.txt:

  job name
  email address
  chain

- input.pdb:

  input PDB file

- input.seq:

  input sequence in fasta format or just sequence


Usage:
======

1. Obtain all dependencies (see above).
2. Modify paths in dependencies/__init__.py


Run by executing this line:
===========================
  python cryptosite.py
