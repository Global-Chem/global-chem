.. _gettingstarted:

Getting Started
===============

This page instructs you on how to get started with global-chem. To begin, make sure you have performed
:ref:`installing global-chem <install>`.

Basic Usage
-----------

The simplest way to start is installing the package and looking through the directory of avaliable variables and retrieve
one of the variables of your choice.

>>> from global_chem import GlobalChem
>>> global_chem = GlobalChem()
>>> print (global_chem.functional_groups_smarts)
>>> {'bromine': '[Br]', 'chlorine': '[Cl]', 'fluorine': '[F]', 'acyl bromide': '[CX3](=[OX1])[Br]', ...}

Other Examples
--------------

Retrieving Amino Acids

>>> from global_chem import GlobalChem
>>> global_chem = GlobalChem()
>>> print (global_chem.amino_acid_side_chains)
>>> {"alanine": "C",  "arginine": "CCCCNC(N)=N", "asparagine": "CCC(N)=O", "aspartic acid": "CC(O)=O", ...}