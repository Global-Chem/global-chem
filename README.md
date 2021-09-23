GlobalChem: A content variable store for Chemistry!
===================================================

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
![Python](https://img.shields.io/badge/python-3.6-blue.svg)
[![Documentation Status](https://readthedocs.org/projects/globalchem/badge/?version=latest)](https://globalchem.readthedocs.io/en/latest/?badge=latest)
![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)
[![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem)
[![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)



Global Chem is a content variable store for cheminformaticians. Variables that we often use in research should be added here. 
The main motivation behind this is to eliminate the development time for when folk build their applications or conduct 
research.

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

Announcements
=============

-   Work has began! April 2020 (added the amino acids)

Using GlobalChem
=====================

GlobalChem, initially, is one class object that have a fixed set of variables. Eventually and hopefully as the package 
grows it will be abstracted out by field i.e protein chemistry, analytical, and so on. 

Installation 
============

GlobalChem is going to be distribute via PyPi and as the content store grows we can expand it to other pieces of software
making it accessible to all regardless of what you use. Alternatively, you could have a glance at the source code and copy/paste
it yourself.

Quick Start
===========

```
from global_chem import GlobalChem
global_chem = GlobalChem()
global_chem.amino_acid_side_chains
global_chem.functional_groups_smiles
```

Variables List
==============
- functional groups (SMILES)
- functional groups (SMARTS)
- amino acids (SMILES)
- common organic solvents (SMILES)
- common organic solvents (SMARTS)
- regex patterns (1 so far Mol2)

Structure of GlobalChem
=======================

Currently, the main subpackages are:

- **global-chem**: globalchem main class. 


Genesis
=======

GlobalChem was created because I noticed I was using the same variable across multiple scripts and figure it would be useful
for folk to have.

- Lead Developer [Suliman sharif](http://sulstice.github.io/)
- Artwork [Elena Chow](http://www.chowelena.com/)

* * * * *

External links
==============


