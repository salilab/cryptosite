Basic usage
***********

Once installed, the CryptoSite protocol can be run by means of a command line
tool (``cryptosite``). Each component of the protocol is also a Python package,
which can be called directly from other Python software
(via ``import cryptosite``).

Overview
========

Running the basic protocol consists of these steps:

#. Create a set of input files specifying the structure to probe, and basic
   CryptoSite parameters.

#. Run ``cryptosite mainer_short`` to calculate structural features and
   prepare input files for
   `AllosMod <https://github.com/salilab/allosmod-lib>`_.

#. `Run AllosMod <https://allosmod.readthedocs.io/en/latest/usage.html#set-up-allosmod-protocol>`_
   given the set of input files.

#. Calculate additional features using the AllosMod output
   (``cryptosite soap_clean``, ``cryptosite pocket_parser``,
    ``cryptosite am_bmi``, ``cryptosite analysis``).

#. Gather together all features into a single file (``cryptosite gatherer``).

#. Predict cryptic binding sites with an SVM using the complete set of features
   as input (``cryptosite predicter``).
