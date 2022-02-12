Global-Chem: Collections of common small molecules and their SMILES/SMARTS to support diverse chemical communities
==================================================================================================================


[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
![Python](https://img.shields.io/badge/python-3.6-blue.svg)
[![Build Status](https://app.travis-ci.com/Sulstice/global-chem.svg?branch=master)](https://app.travis-ci.com/Sulstice/global-chem)
[![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master)
![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)
[![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem)
[![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![downloads](https://img.shields.io/pypi/dm/global-chem)](https://img.shields.io/pypi/dm/global-chem)
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)

Global Chem is an open-source record collection for common and rare chemical lists using IUPAC as input and SMILES/SMARTS as output. As 
mostly needed by myself as I search through chemical infinity. 

I have found these lists written in history to be useful, they come from a variety of different fields but are aggregated 
into the most common format of organic chemists (IUPAC) and the common language of the cheminformatician (SMILES) and for 
pattern matching (SMARTS). 

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

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

To Access Nodes:

```python

from global_chem import GlobalChem

gc = GlobalChem()
gc.check_available_nodes()

>>>
'emerging_perfluoro_alkyls', 'montmorillonite_adsorption', 'common_monomer_repeating_units', 'electrophilic_warheads_for_kinases',
```

Fetch Data from Node:

```python

gc = GlobalChem()
epa = gc.get_node('emerging_perfluoro_alkyls')

epa.get_smiles()
epa.get_smarts()

>>>
{'perfluorohexanoic acid': 'C(=O)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)O'}
```

To Create Your Own Chemical Graph Network
```

from global_chem import GlobalChem

gc = GlobalChem()
gc.check_available_nodes()
gc.initiate_network()

gc.add_node('common_monomer_repeating_units')
gc.add_node('electrophilic_warheads_for_kinases', 'common_monomer_repeating_units')
gc.add_node('common_warhead_covalent_inhibitors', 'common_monomer_repeating_units')
gc.get_parent('common_monomer_repeating_units')    
    

```

Using GlobalChem
=====================

GlobalChem, initially, is one class object with a series of Nodes that are act as objects for any common chemical lists. 
The chemical lists can be accessed as nodes and the user can construct their own node trees for the lists.

An example of it's usage can be found below:

[![asciicast](https://asciinema.org/a/dpbEUdb5SRihzu6uDnRNmVlTv.svg)](https://asciinema.org/a/dpbEUdb5SRihzu6uDnRNmVlTv)

Also since these lists of commonality are stored on github it is easily searchable and tied directly to the paper for 
any bypasser. 

<img width="1440" alt="Screen Shot 2022-01-19 at 9 10 14 AM" src="https://user-images.githubusercontent.com/11812946/150147664-df1149e1-c43b-48b8-946c-4543f39a8bc6.png">

Genesis
=======

GlobalChem was created because I noticed I was using the same variable across multiple scripts and figure it would be useful
for folk to have.

- Lead Developer [Suliman sharif](http://sulstice.github.io/)
- Artwork [Elena Chow](http://www.chowelena.com/)

* * * * *

Citation
========

It's on it's way

