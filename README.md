<h1 align="center">Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities </h1>

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

#### GlobalChem

| Usage     | Build  | Release | SMILES Validation | Licensing | Stats |  
| --------- | ------ | ------- | ----------------- | --------- | ----- |
| [![Documentation](https://img.shields.io/badge/GitBook-Docu-lightblue)](https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/) | [![Test System](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml) | [![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem) | [![saythanks](https://img.shields.io/badge/RDKit-100%25-fg49b4.svg)](https://www.rdkit.org/) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) | [![Downloads](https://pepy.tech/badge/global-chem)](https://pepy.tech/project/global-chem) | 
| [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.chemicalgraphtheory.com) | [![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master) | [![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250) | [![saythanks](https://img.shields.io/badge/PartialSMILES-85.7%25-fg49b4.svg)](https://github.com/baoilleach/partialsmiles) | [![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_shield) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1HaDvAYpaX0_2Eqx1NSl5uJTRkE_0na-L?usp=sharing) | [![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/) | [![Maturity badge - level 2](https://img.shields.io/badge/Maturity-Level%202%20--%20First%20Release-yellowgreen.svg)](https://github.com/tophat/getting-started/blob/master/scorecard.md) | [![saythanks](https://img.shields.io/badge/PySMILES-99.8%25-fg49b4.svg)](https://github.com/pckroon/pysmiles)  | | [![Man Hours](https://img.shields.io/endpoint?url=https%3A%2F%2Fmh.jessemillar.com%2Fhours%3Frepo%3Dhttps%3A%2F%2Fgithub.com%2FSulstice%2Fglobal-chem)](https://jessemillar.com/r/man-hours) |
| ![Python](https://img.shields.io/badge/python-3.6-blue.svg) | [![publish](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml) | | [![saythanks](https://img.shields.io/badge/MolVS-98.5%25-lm89b5.svg)](https://github.com/mcs07/MolVS) | | ![visitor badge](https://visitor-badge.glitch.me/badge?page_id=sulstice.global-chem) |
| | [![pre-commit.ci status](https://results.pre-commit.ci/badge/github/Sulstice/global-chem/master.svg)](https://results.pre-commit.ci/latest/github/Sulstice/global-chem/master) | | [![saythanks](https://img.shields.io/badge/DeepSMILES-99.25%25-lm89b5.svg)](https://github.com/baoilleach/deepsmiles) | |
| | | | [![saythanks](https://img.shields.io/badge/SELFIES-100%25-lm89b5.svg)](https://github.com/aspuru-guzik-group/selfies) | |

#### GlobalChemExtensions

| Release   | Licensing | Stats |
| --------- | --------- | ----- |
| [![PyPI version](https://badge.fury.io/py/global-chem-extensions.svg)](https://badge.fury.io/py/global-chem-extensions) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) | [![Downloads](https://pepy.tech/badge/global-chem-extensions)](https://pepy.tech/project/global-chem-extensions) |

Different extensions demos can be found below, or download them all. Note the PDF extension is the most dependency heavy. 

| Extension              | Demo |
| ---------------------- | ---- |
| forcefields            |      |
| bioinformatics         |      |
| cheminformatics        |      |
| quantum_chemistry      |      |
| development_operations |      |

Installation 
============

GlobalChem is going to be distribute via PyPi and as the tree and it's extensions grows we can expand it to other pieces of software
making it accessible to all regardless of what you use. Alternatively, you could have a glance at the source code and copy/paste
it yourself.

```

pip install global-chem

```

If you want to install the extensions package for extra functionality:


```

pip install global-chem[graphing]
pip install global-chem[forcefields]
pip install global-chem[bioinformatics]
pip install global-chem[cheminformaitcs]
pip install global-chem[quantum_chemistry]
pip install global-chem[development_operations]

```

QuickStart
==========

#### GlobalChem

```

from global_chem import GlobalChem

gc = GlobalChem()

gc.build_global_chem_network(print_output=False, debugger=False)
smiles_list = list(gc.get_node_smiles('pihkal').values())

print (smiles_list)
```

#### GlobalChemExtensions

```

from global_chem_extensions import GlobalChemExtensions

gce = GlobalChemExtensions()

global_chem_molecule = gce.initialize_globalchem_molecule(
    smiles_list[0],
    # stream_file='cgenff.str',
    # frcmod_file='gaff2.frcmod',
)

global_chem_molecule.determine_name()

name = global_chem_molecule.name
attributes = global_chem_molecule.get_attributes()

gce.node_pca_analysis(smiles_list, save_file=False)

```

GlobalChem
==========

### Rules


The Graph Network (GN)s comes with a couple of rules that for now make the software engineering easier on the developer. 

- There must be a root node.
- When Adding a Node every node must be connected. 
- To remove a node it must not have any children. 

The Deep Graph Network (DGN)s comes also with a couple of rules to make the implementation easier:

- There must be a root node of 1 which marks as your "input" node. 
- When adding a layer all nodes will be added to all the previous layers as children. (Folk can use the remove node feature to perform dropouts)

### Knowledge Graph


Just with no dependencies, intialize the class and there you go! All the common and rare groups of the world
at your disposal 

```

gc = GlobalChem()
gc.print_globalchem_network()

                                ┌solvents─common_organic_solvents
             ┌organic_synthesis─└protecting_groups─amino_acid_protecting_groups
             │          ┌polymers─common_monomer_repeating_units
             ├materials─└clay─montmorillonite_adsorption
             │                            ┌privileged_kinase_inhibtors
             │                            ├privileged_scaffolds
             ├proteins─kinases─┌scaffolds─├iupac_blue_book_substituents
             │                 │          └common_r_group_replacements
             │                 └braf─inhibitors
             │              ┌vitamins
             │              ├open_smiles
             ├miscellaneous─├amino_acids
             │              └regex_patterns
global_chem──├environment─emerging_perfluoroalkyls
             │          ┌schedule_one
             │          ├schedule_four
             │          ├schedule_five
             ├narcotics─├pihkal
             │          ├schedule_two
             │          └schedule_three
             ├interstellar_space
             │                    ┌cannabinoids
             │                    │         ┌electrophillic_warheads_for_kinases
             │                    ├warheads─└common_warheads_covalent_inhibitors
             └medicinal_chemistry─│      ┌phase_2_hetereocyclic_rings
                                  └rings─├iupac_blue_book_rings
                                         └rings_in_drugs
                                         

```

### Nodes Contributors


Please follow the node contribution guidelines if you would like to elect your own or someone elses. 

```
'global_chem': Node,
'emerging_perfluoroalkyls': EmergingPerFluoroAlkyls,                      # Asuka Orr & Suliman Sharif
'montmorillonite_adsorption': MontmorilloniteAdsorption,                  # Asuka Orr & Suliman Sharif
'common_monomer_repeating_units': CommonMonomerRepeatingUnits,            # Suliman Sharif
'electrophilic_warheads_for_kinases': ElectrophilicWarheadsForKinases,    # Ruibin Liu & Suliman Sharif
'common_warheads_covalent_inhibitors': CommonWarheadsCovalentInhibitors,  # Shaoqi Zhao & Suliman Sharif
'rings_in_drugs': RingsInDrugs,                                           # Alexander Mackerell & Suliman Sharif
'iupac_blue_book_rings': IUPACBlueBookRings,                              # Suliman Sharif
'phase_2_hetereocyclic_rings': Phase2HetereoCyclicRings,                  # Suliman Sharif
'privileged_scaffolds': PrivilegedScaffolds,                              # Suliman Sharif
'iupac_blue_book': IUPACBlueBook,                                         # Suliman Sharif
'common_r_group_replacements': CommonRGroupReplacements,                  # Sunhwan Jo & Suliman Sharif
'braf_inhibitors': BRAFInhibitors,                                        # Aarion Romany
'privileged_kinase_inhibitors': PrivilegedKinaseInhibitors,               # Suliman Sharif
'common_organic_solvents': CommonOrganicSolvents,                         # Suliman Sharif
'amino_acid_protecting_groups': AminoAcidProtectingGroups,                # Aziza Frank & Suliman Sharif
'schedule_one': ScheduleOne,                                              # Suliman Sharif
'schedule_two': ScheduleTwo,                                              # Suliman Sharif
'schedule_three': ScheduleThree,                                          # Suliman Sharif
'schedule_four': ScheduleFour,                                            # Suliman Sharif
'schedule_five': ScheduleFive,                                            # Suliman Sharif
'interstellar_space': InterstellarSpace,                                  # Suliman Sharif
'vitamins': Vitamins,                                                     # Suliman Sharif
'open_smiles': OpenSmiles,                                                # Suliman Sharif
'amino_acids': AminoAcids,                                                # Suliman Sharif
'pihkal': Pihkal,                                                         # Suliman Sharif
'nickel_ligands': NickelBidendatePhosphineLigands,                        # Suliman Sharif
'cimetidine_and_acyclovir': CimetidineAndAcyclovir,                       # Suliman Sharif
'common_regex_patterns': CommonRegexPatterns,                             # Chris Burke & Suliman Sharif
```

<details><summary><h3>Node List</h1><br/></summary>

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
| Excipients Cimetidine & Acyclovir   | 14           | Vaithianathan, Soundarya, et al. “Effect of Common Excipients on the Oral Drug Absorption of Biopharmaceutics Classification System Class 3 Drugs Cimetidine and Acyclovir.” Journal of Pharmaceutical Sciences, vol. 105, no. 2, Feb. 2016, pp. 996–1005. PubMed, https://doi.org/10.1002/jps.24643.                |
| Nickel Bidendate Phosphine Ligands  | N/A          | Clevenger, Andrew L., et al. “Trends in the Usage of Bidentate Phosphines as Ligands in Nickel Catalysis.” Chemical Reviews, vol. 120, no. 13, July 2020, pp. 6124–96. DOI.org (Crossref), https://doi.org/10.1021/acs.chemrev.9b00682.                                                                              |
| Common Regex Patterns               | 1            |                                                                                                     

</details>                            

GlobalChemExtensions
====================

A Variety of Tools are available for you to browse and analyze data and with the full list of different applications can be found in the google colab demo or the Gitbook documentation. A demonstration of the data visualization extensions designed with plotly and bokeh are displayed below:

<p align="center">
  <img width="800" height="600" src="https://raw.githubusercontent.com/Sulstice/global-chem/master/images/figures/figure_10.png">
</p>

<details><summary><h3>Extension List</h1><br/></summary>

| Extension                       | Description                                                                                                             | Feature Import   |
|---------------------------------|-------------------------------------------------------------------------------------------------------------------------|------------------|
| GlobalChem Chemical Entities    | GlobalChem has internal Molecule objects with all common attributes associated and conversion to SMILES                 | extensions       |
| GlobalChem Biological Entities  | GlobalChem has internal DNA/RNA/Protein/Molecule objects with all common attributes associated and conversion to SMILES | bioinformatics   |
| ForceField Molecules            | GlobalChem can parse, manipulate, and write CGenFF and GaFF2 files as objects                                           | extensions       |
| PDF Generation and Parsing      | GlobalChem can generate SMILES to PDF and convert the PDF to SMILES                                                     | pdf              |
| SMILES Validation               | GlobalChem has connection to PySMILES, DeepSMILES, PartialSmiles, SELFIES, MolVS for validation of SMILES sets          | validation       |
| SMILES Protonation States       | GlobalChem can take a set of compounds and predict the protonation states of a SMILES string over a range of pH         | extensions       |
| Open Source Database Monitoring | GlobalChem uses Uptime-Cheminformatics to Keep Track of Open Source Chemical Data                                       | extensions       |
| Networkx Software Adapter       | GlobalChem Network can be converted into NetworkX Graph Objects                                                         | extensions       |
| SMARTS Pattern Validation       | GlobalChem uses the MiniFrag Database to test SMARTS strings accuracy for functional group selection                    | extensions       |
| Principal Component Analysis    | GlobalChem can readily interpret SMILES, fingerprint, cluster and apply PCA analysis user can tweak parameters          | machine_learning |
| Drug Design Filters             | GlobalChem can filter compounds based on Common Drug Design Filtering Rules                                            | extensions       |
| Deep Layer Scatter Analysis     | To visualize relations between sets of molecules, GlobalChem offers a parallel coordinate diagram generation            | machine_learning | 
| Sunbursting Radial Analysis     | GlobalChem offers a sunbursting mechanism to allow uses to observe how sets of compounds relate to the common set      | machine_learning |
| Graphing Templates              | GlobalChem offers graphing templates to aid in faster data analysis, currently the only offer is Plotly               | machine_learning |
| CGenFF Dissimilarity Score      | GlobalChem can offer the difference between two molecules based on their Atom Types                                     | extensions       |
| OneHot Encoding                 | GlobalChem has it's own one hot encoder and decoder based on the common lists for Machine Learning                      | extensions       |
| SMARTS Pattern Identifier       | GlobalChem connects to the SMARTS Plus and can offer visualization into different SMARTS components                     | extensions       |

</details>

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

* * * * *

Citation
========

It's on it's way


Licensing
=========

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)

