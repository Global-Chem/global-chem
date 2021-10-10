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



Global Chem is a variable store for common and rare chemical lists using IUPAC as input and SMILES/SMARTS as output. As 
mostly needed by myself as I search through chemical infinity. 

I have found these lists written history to be useful, they come from a variety of different fields but are aggregated 
into the most common format of organic chemists (IUPAC) and the common language of the cheminformatician (SMILES) and for 
pattern matching (SMARTS).

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>



Using GlobalChem
=====================

GlobalChem, initially, is one class object that have a fixed set of variables. Eventually and hopefully as the package 
grows it will be abstracted out by field i.e protein chemistry, analytical, and so on. 

Installation 
============

GlobalChem is going to be distribute via PyPi and as the content store grows we can expand it to other pieces of software
making it accessible to all regardless of what you use. Alternatively, you could have a glance at the source code and copy/paste
it yourself.

```

pip install global-chem

```
Quick Start
===========

Just with no dependencies, intialize the class and there you go! All the common and rare groups of the world
at your disposal 

```

from global_chem import GlobalChem

global_chem = GlobalChem()


```

Variables List
==============

| Chemical List                       | Languages                    | Variables                                                                                                                  | References                                                                                                                                                                                                                                     |
|-------------------------------------|------------------------------|----------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | amino_acid_side_smiles, amino_acid_side_smarts                                                                             | Common Knowledge                                                                                                                                                                                                                               |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | common_organic_solvents_smiles, common_organic_solvents_smarts                                                             | Fulmer, Gregory R., et al. “NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist.”Organometallics , vol. 29, no. 9, May 2010, pp. 2176–79.   |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | open_smiles_functional_groups_smiles, open_smiles_functional_groups_smarts                                                 | OpenSMILES Home Page. http://opensmiles.org/.                                                                                                                                                                                                  |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | iupac_blue_book_radical_smiles, iupac_blue_book_radical_smarts, iupac_blue_book_rings_smiles, iupac_blue_book_rings_smarts | Chemical Rubber Company. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data . Edited by David R. Lide, 85. ed, CRC Press, 2004.                                                                       |
| Common Regex Patterns               | Mol2                         | common_regex_patterns                                                                                                      |                                                                                                                                                                                                                                                |


Genesis
=======

GlobalChem was created because I noticed I was using the same variable across multiple scripts and figure it would be useful
for folk to have.

- Lead Developer [Suliman sharif](http://sulstice.github.io/)
- Artwork [Elena Chow](http://www.chowelena.com/)

* * * * *

Citation
========

If you find my personal dictionary is useful please feel free to cite me under:

```

Sul, Elena Y. Chow: "Sul's Dictionary" 

```
External links
==============


