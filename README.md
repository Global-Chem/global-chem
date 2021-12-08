GlobalChem: A content variable store for Chemistry!
===================================================

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
![Python](https://img.shields.io/badge/python-3.6-blue.svg)
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

compounds = global_chem.rings_in_drugs_smiles

```

Variables List
==============

| Chemical List                       | Languages                    | Variables                                                                                                                  | # of Entries | References                                                                                                                                                                                                                                                                                                           |
|-------------------------------------|------------------------------|----------------------------------------------------------------------------------------------------------------------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | amino_acid_side_smiles, amino_acid_side_smarts                                                                             | 20           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Essential Vitamins                  | Preferred Name/SMILES/SMARTS | vitamins_smiles, vitamins_smarts                                                                                           | 13           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | common_organic_solvents_smiles, common_organic_solvents_smarts                                                             | 42           | Fulmer, Gregory R., et al. “NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist.”Organometallics, vol. 29, no. 9, May 2010, pp. 2176–79.                                                                          |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | open_smiles_functional_groups_smiles, open_smiles_functional_groups_smarts                                                 | 94           | OpenSMILES Home Page. http://opensmiles.org/.                                                                                                                                                                                                                                                                        |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | iupac_blue_book_radical_smiles, iupac_blue_book_radical_smarts, iupac_blue_book_rings_smiles, iupac_blue_book_rings_smarts | 333          | Chemical Rubber Company. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edited by David R. Lide, 85. ed, CRC Press, 2004.                                                                                                                                               |
| Rings in Drugs                      | IUPAC/SMILES/SMARTS          | rings_in_drugs_smiles, rings_in_drugs_smarts                                                                               | 92           | Taylor, Richard D., et al. “Rings in Drugs.” Journal of Medicinal Chemistry, vol. 57, no. 14, July 2014, pp. 5845–59. ACS Publications, https://doi.org/10.1021/jm4017625.                                                                                                                                           |
| Phase 2 Hetereocyclic Rings         | IUPAC/SMILES/SMARTS          | phase_2_hetereocyclic_rings_smiles, phase_2_hetereocyclic_rings_smarts                                                     | 19           | Broughton, Howard B., and Ian A. Watson. “Selection of Heterocycles for Drug Design.” Journal of Molecular Graphics & Modelling, vol. 23, no. 1, Sept. 2004, pp. 51–58. PubMed, https://doi.org/10.1016/j.jmgm.2004.03.016.                                                                                          |
| Privileged Scaffolds                | IUPAC/SMILES/SMARTS          | privileged_scaffolds_smiles, privileged_scaffolds_smarts                                                                   | 47           | Welsch, Matthew E., et al. “Privileged Scaffolds for Library Design and Drug Discovery.” Current Opinion in Chemical Biology , vol. 14, no. 3, June 2010, pp. 347–61.PubMed, https://doi.org/10.1016/j.cbpa.2010.02.018.                                                                                             |
| Common Warheads                     | IUPAC/SMILES/SMARTS          | common_warhead_smiles, common_warhead_smarts                                                                               | 29           | Gehringer, Matthias, and Stefan A. Laufer. “Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology.”Journal of Medicinal Chemistry , vol. 62, no. 12, June 2019, pp. 5673–724. ACS Publications, https://doi.org/10.1021/acs.jmedchem.8b01153. |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | common_polymer_repeating_units_smiles, common_polymer_repeating_units_smarts                                               | 78           | Hiorns, R. C., et al. “A brief guide to polymer nomenclature (IUPAC Technical Report).”Pure and Applied Chemistry , vol. 84, no. 10, Oct. 2012, pp. 2167–69., https://doi.org/10.1351/PAC-REP-12-03-05.                                                                                                              |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | r_groups_replacements_smiles, r_groups_replacements_smarts                                                                 | 499          | Takeuchi, Kosuke, et al. “R-Group Replacement Database for Medicinal Chemistry.”   Future Science OA , vol. 7, no. 8, Sept. 2021, p. FSO742.   future-science.com (Atypon) , https://doi.org/10.2144/fsoa-2021-0062.                                                                                                 |
| Common Regex Patterns               | Mol2                         | common_regex_patterns                                                                                                      | 1            |                                                                                                                                                                                                                                                                                                                      |

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


