Installation
************

In the Sali Lab
===============

If you are working in the Sali lab, you don't need to build and install
CryptoSite - it is already set up for you as a module. Just run
``module load cryptosite`` to load it.

Dependencies
============

All dependencies listed below are expected to be found in standard
system paths. This may require setting ``PYTHONPATH``, ``PATH`` and/or
``LD_LIBRARY_PATH`` environment variables. Note that Linux is the only platform
for which all these dependencies are available, and so is the only platform
on which CryptoSite currently functions.

* `Python <https://www.python.org>`_ 2.6 or later (Python 3 should also be OK).

* `MODELLER <https://salilab.org/modeller/>`_ plus the
  `SOAP-Protein library <https://salilab.org/SOAP/>`_.
  The ``soap_protein_od.hdf5`` file needs to be placed into MODELLER's
  ``modlib`` directory.
 
* `MUSCLE <http://www.drive5.com/muscle/>`_.

* `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_. It is expected that the
  :command:`mkdssp` binary is in the ``PATH``.

* `fpocket <http://fpocket.sourceforge.net/>`_ (version 2).

* `PatchDock <http://bioinfo3d.cs.tau.ac.il/PatchDock/>`_.

* `IMP <https://integrativemodeling.org/>`_.

* `NCBI BLAST+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_
  plus a local copy of the `UniProt database <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/>`_ for it to search against.
  To make this local database, ``gunzip`` the ``uniprot_trembl.fasta.gz``
  and ``uniprot_sprot.fasta.gz`` files available from the UniProt website
  and combine them into a single ``uniprot`` text file.

* `USEARCH <http://drive5.com/usearch/download.html>`_ version 8.1 or later.

* `Biopython <http://biopython.org/>`_.

* `NumPy <http://www.numpy.org/>`_ and `SciPy <https://scipy.org/scipylib/>`_.

* `scikit-learn <http://scikit-learn.org/>`_. Note that precisely version 0.12
  is needed - other versions won't work.

* `AllosMod <https://github.com/salilab/allosmod-lib>`_ is needed to run part
  of the protocol.

* `nose <https://nose.readthedocs.io/en/latest/>`_ is also needed to run the
  test suite (recommended but not essential).

In the Sali lab, running 
``module load modeller muscle dssp fpocket patch_dock imp blast+ usearch``
will get all of these dependencies.

Building
========

Use ``make PYTHON=python3`` or ``make PYTHON=python2`` to build the library
(depending on which version of Python you want to use).
Use ``make test`` to test the library, and ``make install`` to install it.
In most cases you will need to tell ``make`` where to install (if running on
a Linux cluster, CryptoSite will need to be installed on a network-accessible
filesystem) and where your local copy of UniProt is, with something like
``make PREFIX=/shared/cryptosite UNIPROT=/database/uniprot install``. See
``Makefile.include`` for all make variables that can be configured.
