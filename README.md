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

Using GlobalChem
=====================

GlobalChem, initially, is one class object with a series of Nodes that are act as objects for any common chemical lists. 
The chemical lists can be accessed as nodes and the user can construct their own node trees for the lists.

An example of it's usage can be found below:


Also since these lists of commonality are stored on github it is easily searchable and tied directly to the paper for 
any bypasser. 

<img width="1440" alt="Screen Shot 2022-01-19 at 9 10 14 AM" src="https://user-images.githubusercontent.com/11812946/150147664-df1149e1-c43b-48b8-946c-4543f39a8bc6.png">

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

Variables List
==============

| Chemical List                       | Languages                    | Variables                                            | # of Entries | References                                                                                                                                                                                                                                                                                                           |
|-------------------------------------|------------------------------|------------------------------------------------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | get_amino_acids()                                    | 20           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Essential Vitamins                  | Preferred Name/SMILES/SMARTS | get_essential_vitamins()                             | 13           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | get_common_organic_solvents()                        | 42           | Fulmer, Gregory R., et al. “NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist.”Organometallics, vol. 29, no. 9, May 2010, pp. 2176–79.                                                                          |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | get_open_smiles_functional_groups()                  | 94           | OpenSMILES Home Page. http://opensmiles.org/.                                                                                                                                                                                                                                                                        |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | get_iupac_blue_book_common_functional_groups()       | 333          | Chemical Rubber Company. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edited by David R. Lide, 85. ed, CRC Press, 2004.                                                                                                                                               |
| Rings in Drugs                      | IUPAC/SMILES/SMARTS          | get_rings_in_drugs()                                 | 92           | Taylor, Richard D., et al. “Rings in Drugs.” Journal of Medicinal Chemistry, vol. 57, no. 14, July 2014, pp. 5845–59. ACS Publications, https://doi.org/10.1021/jm4017625.                                                                                                                                           |
| Phase 2 Hetereocyclic Rings         | IUPAC/SMILES/SMARTS          | get_common_heterocyclic_rings_phase_2()              | 19           | Broughton, Howard B., and Ian A. Watson. “Selection of Heterocycles for Drug Design.” Journal of Molecular Graphics & Modelling, vol. 23, no. 1, Sept. 2004, pp. 51–58. PubMed, https://doi.org/10.1016/j.jmgm.2004.03.016.                                                                                          |
| Privileged Scaffolds                | IUPAC/SMILES/SMARTS          | get_common_privileged_scaffolds()                    | 47           | Welsch, Matthew E., et al. “Privileged Scaffolds for Library Design and Drug Discovery.” Current Opinion in Chemical Biology , vol. 14, no. 3, June 2010, pp. 347–61.PubMed, https://doi.org/10.1016/j.cbpa.2010.02.018.                                                                                             |
| Common Warheads                     | IUPAC/SMILES/SMARTS          | get_common_warhead_covalent_inhibitors()             | 29           | Gehringer, Matthias, and Stefan A. Laufer. “Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology.”Journal of Medicinal Chemistry , vol. 62, no. 12, June 2019, pp. 5673–724. ACS Publications, https://doi.org/10.1021/acs.jmedchem.8b01153. |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | get_common_polymer_repeating_units()                 | 78           | Hiorns, R. C., et al. “A brief guide to polymer nomenclature (IUPAC Technical Report).”Pure and Applied Chemistry , vol. 84, no. 10, Oct. 2012, pp. 2167–69., https://doi.org/10.1351/PAC-REP-12-03-05.                                                                                                              |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | get_commonly_used_r_group_replacements()             | 499          | Takeuchi, Kosuke, et al. “R-Group Replacement Database for Medicinal Chemistry.”   Future Science OA , vol. 7, no. 8, Sept. 2021, p. FSO742.   future-science.com (Atypon) , https://doi.org/10.2144/fsoa-2021-0062.                                                                                                 |
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | get_common_electrophilic_warheads_for_kinases()      | 24           | Petri, László, et al. “An Electrophilic Warhead Library for Mapping the Reactivity and Accessibility of Tractable Cysteines in Protein Kinases.” European Journal of Medicinal Chemistry, vol. 207, Dec. 2020, p. 112836. PubMed, https://doi.org/10.1016/j.ejmech.2020.112836.                                      |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | get_privileged_scaffolds_for_kinase_inhibitors()     | 29           | Hu, Huabin, et al. “Systematic Comparison of Competitive and Allosteric Kinase Inhibitors Reveals Common Structural Characteristics.” European Journal of Medicinal Chemistry, vol. 214, Mar. 2021, p. 113206. ScienceDirect, https://doi.org/10.1016/j.ejmech.2021.113206.                                          |
| BRaf Inhibitors                     | IUPAC/SMILES/SMARTS          | get_braf_kinase_inhibitors_for_cancer()              | 54           | Agianian, Bogos, and Evripidis Gavathiotis. “Current Insights of BRAF Inhibitors in Cancer.” Journal of Medicinal Chemistry, vol. 61, no. 14, July 2018, pp. 5775–93. ACS Publications, https://doi.org/10.1021/acs.jmedchem.7b01306.                                                                                |
| Common Amino Acid Protecting Groups | IUPAC/ACRONYM/SMILES/SMARTS  | get_common_amino_acid_protecting_groups()            | 346          | Isidro-Llobet, Albert, et al. “Amino Acid-Protecting Groups.” Chemical Reviews, vol. 109, no. 6, June 2009, pp. 2455–504. DOI.org (Crossref), https://doi.org/10.1021/cr800323s.                                                                                                                                     |
| Emerging Perfluoroalkyls            | IUPAC/SMILES/SMARTS          | get_polyfluoroalkyl_substances()                     | 27           | Pelch, Katherine E., et al. “PFAS Health Effects Database: Protocol for a Systematic Evidence Map.” Environment International, vol. 130, Sept. 2019, p. 104851. ScienceDirect, https://doi.org/10.1016/j.envint.2019.05.045.                                                                                         |
| Chemicals For Clay Adsorption       | IUPAC/SMILES/SMARTS          | get_chemical_adsorption_on_montmorillonite_clays()   | 33           | Orr, Asuka A., et al. “Combining Experimental Isotherms, Minimalistic Simulations, and a Model to Understand and Predict Chemical Adsorption onto Montmorillonite Clays.” ACS Omega, vol. 6, no. 22, June 2021, pp. 14090–103. PubMed, https://doi.org/10.1021/acsomega.1c00481.                                     |
| Cannabinoids                        | IUPAC/SMILES/SMARTS          | get_cannabinoid_smiles(), get_cannabinoid_smarts()   | 63           | Turner, Carlton E., et al. “Constituents of Cannabis Sativa L. XVII. A Review of the Natural Constituents.” Journal of Natural Products, vol. 43, no. 2, Mar. 1980, pp. 169–234. ACS Publications, https://doi.org/10.1021/np50008a001.                                                                              |
| Schedule 1 United States Narcotics  | Preferred Name/SMILES/SMARTS | get_schedule_one()                                   | 240          | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 2 United States Narcotics  | Preferred Name/SMILES/SMARTS | get_schedule_two()                                   | 60           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 3 United States Narcotics  | Preferred Name/SMILES/SMARTS | get_schedule_three()                                 | 22           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 4 United States Narcotics  | Preferred Name/SMILES/SMARTS | get_schedule_four()                                  | 77           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 5 United States Narcotics  | Preferred Name/SMILES/SMARTS | get_schedule_five()                                  | 8            | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Common Regex Patterns               | Mol2                         | common_regex_patterns                                | 1            |                                                                                                                                                                                                                                                                                                                      |

Data Collection
===============

References and associatied compound lists are selected based on the interests of the scientific contributors.  This should include consideration of relevance to the scientific community. 
The SMILES strings may be abstracted in a variety of methods:

-  For simple molecules one representation of the SMILES can be directly translated using visual 
inspection. This is typically appropriate for compounds at the beginning of a reported list that contain the most common denominator rings. 

- For complex molecules the image can be redrawn in the free version of ChemDraw and then translated into SMILES. 

- For sources where the SMILES are written and the IUPAC is not known the SMILES are translated into ChemDraw and the name retrieved. 
Note that some of the names may be modified based on human inspection in favor of preferred names. 

- For polymer papers, the site points were omitted from the name and some of the nomenclature adjusted for preferred names
over traditional. For example: 'yl' to mark site points for polymer connections was removed in favor of reduced english complexity. 

- In the case of radicals, some SMILES were adjusted to remove the radical chemical feature as they serve as connection points. However in some cases the radical component was maintained, especially in the case of IUPAC blue book common substituents.

- SMARTS strings were adapted from the SMILES using RDKit [@Landrum:2019-5]


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

