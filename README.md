Global-Chem: A Chemical Graph Network of common small molecules and their SMILES/SMARTS to support diverse chemical communities
===============================================================================================================================


[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
![Python](https://img.shields.io/badge/python-3.6-blue.svg)
[![Build Status](https://app.travis-ci.com/Sulstice/global-chem.svg?branch=master)](https://app.travis-ci.com/Sulstice/global-chem)
[![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master)
![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)
[![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem)
[![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![Downloads](https://pepy.tech/badge/global-chem)](https://pepy.tech/project/global-chem)
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_shield)

Global Chem is an open-source graph database and api for common and rare chemical lists using IUPAC as input and SMILES/SMARTS as output. As 
mostly needed by myself as I search through chemical infinity.

I have found these lists written in history to be useful, they come from a variety of different fields but are aggregated 
into the most common format of organic chemists (IUPAC) and the common language of the cheminformatician (SMILES) and for 
pattern matching (SMARTS). 

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

Overview
========

|      |      |
| ---- | ---- |
| Link to Documentation | [![Documentation](https://img.shields.io/badge/GitBook-Docu-lightblue)](https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/) |
| Link to Website | [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.chemicalgraphtheory.com)

Installation 
============

GlobalChem is going to be distribute via PyPi and as the content store grows we can expand it to other pieces of software
making it accessible to all regardless of what you use. Alternatively, you could have a glance at the source code and copy/paste
it yourself.

```

pip install global-chem

```

If you want to install the extensions package for extra functionality but with dependencies: https://github.com/Sulstice/global-chem-extensions

```

pip install global-chem-extensions

```

Rules
=====

The Graph Network (GN)s comes with a couple of rules that for now make the software engineering easier on the developer. 

1.) There must be a root node.
2.) When Adding a Node every node must be connected. 
3.) To remove a node it must not have any children. 

The Deep Graph Network (DGN)s comes also with a couple of rules to make the implementation easier:

1.) There must be a root node of 1 which marks as your "input" node. 
2.) When adding a layer all nodes will be added to all the previous layers as children. (Folk can use the remove node feature to perform dropouts)

Quick Start
===========

Just with no dependencies, intialize the class and there you go! All the common and rare groups of the world
at your disposal 

Adding Your Own Chemical List
=============================

If you would like to add your paper to the chemical graph network then please "File an Issue" with your chemical list 
and perhaps a suggestion of where to add it or you can leave for up to us to decide. The format of the chemical list can
be something like this:

```

smiles = {
   '3,5-dimethoxyphenylisoproxycarbonyl': 'COC1=CC(C(C)(OC=O)C)=CC(OC)=C1',
   '2-(4-biphenyl)isopropoxycarbonyl': 'CC(C)(OC=O)C(C=C1)=CC=C1C2=CC=CC=C2',
   '2-nitrophenylsulfenyl': 'SC1=CC=CC=C1[N+]([O-])=O',
   'boc': 'O=COC(C)(C)C',
}  

```

Nodes List
==============

| Chemical List                       | # of Entries | References                                                                                                                                                                                                                                                                                                           |
|-------------------------------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Amino Acids                         | 20           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Essential Vitamins                  | 13           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Common Organic Solvents             | 42           | Fulmer, Gregory R., et al. “NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist.”Organometallics, vol. 29, no. 9, May 2010, pp. 2176–79.                                                                          |
| Open Smiles                         | 94           | OpenSMILES Home Page. http://opensmiles.org/.                                                                                                                                                                                                                                                                        |
| IUPAC Blue Book (CRC Handbook) 2003 | 333          | Chemical Rubber Company. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edited by David R. Lide, 85. ed, CRC Press, 2004.                                                                                                                                               |
| Rings in Drugs                      | 92           | Taylor, Richard D., et al. “Rings in Drugs.” Journal of Medicinal Chemistry, vol. 57, no. 14, July 2014, pp. 5845–59. ACS Publications, https://doi.org/10.1021/jm4017625.                                                                                                                                           |
| Phase 2 Hetereocyclic Rings         | 19           | Broughton, Howard B., and Ian A. Watson. “Selection of Heterocycles for Drug Design.” Journal of Molecular Graphics & Modelling, vol. 23, no. 1, Sept. 2004, pp. 51–58. PubMed, https://doi.org/10.1016/j.jmgm.2004.03.016.                                                                                          |
| Privileged Scaffolds                | 47           | Welsch, Matthew E., et al. “Privileged Scaffolds for Library Design and Drug Discovery.” Current Opinion in Chemical Biology , vol. 14, no. 3, June 2010, pp. 347–61.PubMed, https://doi.org/10.1016/j.cbpa.2010.02.018.                                                                                             |
| Common Warheads                     | 29           | Gehringer, Matthias, and Stefan A. Laufer. “Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology.”Journal of Medicinal Chemistry , vol. 62, no. 12, June 2019, pp. 5673–724. ACS Publications, https://doi.org/10.1021/acs.jmedchem.8b01153. |
| Common Polymer Repeating Units      | 78           | Hiorns, R. C., et al. “A brief guide to polymer nomenclature (IUPAC Technical Report).”Pure and Applied Chemistry , vol. 84, no. 10, Oct. 2012, pp. 2167–69., https://doi.org/10.1351/PAC-REP-12-03-05.                                                                                                              |
| Common R Group Replacements         | 499          | Takeuchi, Kosuke, et al. “R-Group Replacement Database for Medicinal Chemistry.”   Future Science OA , vol. 7, no. 8, Sept. 2021, p. FSO742.   future-science.com (Atypon) , https://doi.org/10.2144/fsoa-2021-0062.                                                                                                 |
| Electrophillic Warheads for Kinases | 24           | Petri, László, et al. “An Electrophilic Warhead Library for Mapping the Reactivity and Accessibility of Tractable Cysteines in Protein Kinases.” European Journal of Medicinal Chemistry, vol. 207, Dec. 2020, p. 112836. PubMed, https://doi.org/10.1016/j.ejmech.2020.112836.                                      |
| Privileged Scaffolds for Kinases    | 29           | Hu, Huabin, et al. “Systematic Comparison of Competitive and Allosteric Kinase Inhibitors Reveals Common Structural Characteristics.” European Journal of Medicinal Chemistry, vol. 214, Mar. 2021, p. 113206. ScienceDirect, https://doi.org/10.1016/j.ejmech.2021.113206.                                          |
| BRaf Inhibitors                     | 54           | Agianian, Bogos, and Evripidis Gavathiotis. “Current Insights of BRAF Inhibitors in Cancer.” Journal of Medicinal Chemistry, vol. 61, no. 14, July 2018, pp. 5775–93. ACS Publications, https://doi.org/10.1021/acs.jmedchem.7b01306.                                                                                |
| Common Amino Acid Protecting Groups | 346          | Isidro-Llobet, Albert, et al. “Amino Acid-Protecting Groups.” Chemical Reviews, vol. 109, no. 6, June 2009, pp. 2455–504. DOI.org (Crossref), https://doi.org/10.1021/cr800323s.                                                                                                                                     |
| Emerging Perfluoroalkyls            | 27           | Pelch, Katherine E., et al. “PFAS Health Effects Database: Protocol for a Systematic Evidence Map.” Environment International, vol. 130, Sept. 2019, p. 104851. ScienceDirect, https://doi.org/10.1016/j.envint.2019.05.045.                                                                                         |
| Chemicals For Clay Adsorption       | 33           | Orr, Asuka A., et al. “Combining Experimental Isotherms, Minimalistic Simulations, and a Model to Understand and Predict Chemical Adsorption onto Montmorillonite Clays.” ACS Omega, vol. 6, no. 22, June 2021, pp. 14090–103. PubMed, https://doi.org/10.1021/acsomega.1c00481.                                     |
| Cannabinoids                        | 63           | Turner, Carlton E., et al. “Constituents of Cannabis Sativa L. XVII. A Review of the Natural Constituents.” Journal of Natural Products, vol. 43, no. 2, Mar. 1980, pp. 169–234. ACS Publications, https://doi.org/10.1021/np50008a001.                                                                              |
| Schedule 1 United States Narcotics  | 240          | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 2 United States Narcotics  | 60           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 3 United States Narcotics  | 22           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 4 United States Narcotics  | 77           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 5 United States Narcotics  | 8            | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Pihkal                              | 179          | Shulgin, Alexander T., and Ann Shulgin. Pihkal: A Chemical Love Story. 1. ed., 8. print, Transform, 2010.                                                                                                                                                                                                            |
| Common Regex Patterns               | 1            |                                                                                                                                                                                                                                                                                                                      |

Open Source Software Compliance
===============================

`GlobalChem` follows the same principles outlined in part 11 of Title 21 of the Code of Federal Regulations; Electronic Records,
Electronic Signatures (21 CFR Part 11) guidance documentation. Since there are no formal guidelines for how open source software should be handled, we
attempt at completing requirements. The FDA considers part 11 to be applicable to the following criteria of electronic records and how 
`GlobalChem` accomplishes each component:

- **Plausabilitiy:** `GlobalChem` was built on data that was abstracted from books and papers using reading and redrawing. It adds a component of
IUPAC/SMILES/SMARTS strings to store it electronically which give it's data it's unique component. The records are open sourced
and appropiately version controlled by maintainers of the repository and open source community feedback. 
`GlobalChem`'s purposes are still unknown as it enters open source deployment. We have built extended functions that live in
a seperate package `GlobalChemExtensions` that do depend on `GlobalChem`. Since each version is packaged appropiately, if 
reliance on a version is a need then it's software is available on `Github` and `PyPi`. A Standard Operating Procedure (SOP)
can be filed submitted from the extensions utility documentation maintained on `Gitbook`

- **Validation:** `GlobalChem` follows Good Automated Manufacturing Practice (GAMP) Category 3 which is "software that is used as installed"
and potentially "configurable". `GlobalChem` testing comes from within, the documentation serves as the ultimate test
for functionality because that is what the users will test the most since we rely on open source. A continous integration (CI)
system is also built concomitantly to serve as basic functionality testing of the `GlobalChem` graph network. The Data stored
is maintained by experts in the field but subject to change based on community feedback if an error is found. 

- **Audit Trail:** `GlobalChem` is version controlled with `Git` and hosted on Microsoft's platform `Github`. `GlobalChem` follows a semantic
versioning control of the schema `X1.X2.X3`: `X1` marks formal stable releases with tests and docuementation and mean
big refactoring to the software or in functionality, `X2` means a new feature is added with or without tests and documentation but 
iterates as so. `X3` means a "hot" fix (something that is a an easy bug), small feature or additional parameter to add to a function
, or iteration to the data. 

- **Legacy Systems:** `GlobalChem` has been operational for nearly 2 years since it's first release with version `0.3.0` in May 2020. `GlobalChem`
was built with a full trail in the open source community with each version catalogued and visibility to all. This satisfies 
the rules outlines for determining a legacy system. We use community feedback provided from social media platforms (Twitter, Github, LinkedIn)
as documented evidence and justification that `GlobalChem` is fit for it's intended use of cheminformatics.

- **Copies of Records:** `GlobalChem` has records stored on `Github` for the software that can be exported to a variety of formats as provided by
Microsoft. For documentation, it is hosted on `Gitbook` and versioning controlled in accordance to the software. Each "book"
can be exported into Portable Data Format (PDF) appropiate for FDA submission.

- **Record Retention:** `GlobalChem` has a record of the documentation versioned controlled to a unique id (UUID) that serves as it's identifier
for each iteration stored on `Gitbook`. Each version is stored as markdown files and be converted to PDF, if needed. 

`GlobalChem` has a Mozilla Public License version 2.0. `GlobalChem` allows you to use the software in your larger work and 
extend it with modifications if you wish. The contingency is that if you install `GlobalChem` and release new software 
then you must follow the same principles installed in our license for the open source community. 

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

- SMARTS strings were adapted from the SMILES using RDKit (4)

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


## License
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)
