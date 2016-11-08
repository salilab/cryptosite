Basic usage
***********

Once installed, the CryptoSite protocol can be run by means of a command line
tool (``cryptosite``). Each component of the protocol is also a Python package,
which can be called directly from other Python software
(via ``import cryptosite``).

Overview
========

Running the basic protocol consists of three steps:

#. Create a set of input files specifying the structure to probe, and basic
   CryptoSite parameters.

#. Run ``cryptosite mainer_short`` to calculate structural features and
   prepare input files for
   `AllosMod <https://github.com/salilab/allosmod-lib>`_.

#. `Run AllosMod <https://allosmod.readthedocs.io/en/latest/usage.html#set-up-allosmod-protocol>`_
   given the set of input files.

#. Calculate additional features using the AllosMod output.

#. Predict cryptic binding sites with an SVM using the complete set of features
   as input (``cryptosite predicter``).