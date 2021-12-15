---
title: 'Global-Chem: A record collection of common small molecules and their SMILES/SMARTS in different chemical communities'
tags:
  - Python
  - Cheminformatics
authors:
  - name: Suliman Sharif,  Elena Yi Chow, Asuka Orr, Aarion Romany, Aziza Frank, Shaoqi Zhan, Ruibin Liu, Sunhwan Jo, Alexander D. MacKerell Jr. 
    orcid: 0000-0002-1342-9258
    affiliation: 1
affiliations:
 - name: University of Maryland, School of Pharmacy
   index: 1
date: 12/08/2021
bibliography: paper.bib
---

# Introduction

The in silico chemical universe is expanding rapidly as open access titan databases (Enamine Database (20 Billion) [@Gorgulla:2020-4],
Zinc Database (2 Billion) [Irwin:2020-12], PubMed Database (68 Million) [Roberts:2001-2]) and cheminformatic tools
to process, manipulate, and derive new compound structures are established. This left us with a chemical data big bang
with ultra-large datasets and an ambiguous classification system in an attempt to organize the data. Previously, partial
organizational attempts were made on PubMed filling chemical data linkages for computational toxicology called Actor for a specific
refactored and refined effort [Judson:2019-9]. For the EnamineDB, a scaffold to biological activity was designed to target 
Toll-Like Receptors in an object-oriented fashion [Perez-Regidor:2016-9]. These organizational methods are difficult
to reproduce as well as can be difficult to implement given the amount of data. When applying these papers they don't provide
so much use to the common developer. So what do we do?

To organize the data we need to revert back to the idea of simplistic communication and distribution. Humans use symbols and drawings to communicate, a set of symbols and their combinations
are called a language. Different languages can be employed to carry different features and mean different things to a variety of communities as they infer meaning. 
IUPAC was a written language that predates even drawing atoms as a method of communication between chemists [Cooke-Fox:1989-5]; 
other chemical sub-communities also adopted the language and applied to their field to different dialects i.e polymer chemistry, organo-metallic chemistry.
In the recent years, SMILES [Weininger:1988-5] is becoming a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D or 3D geometry with ease.
Due to it's "first to market" scientific chemical language IUPAC is the legacy language that is a lexical key to unlocking informational wealth about a chemical pattern or group. But there are problems with the language due to it's length in describing bigger molecules. IUPAC names in organic chemisty papers can extend pages with no real value. To compact information, chemists just released the drawings but that can be hard to store precisely. Algorithms
are being designed to abstract and interpolate skeletal patterns and languages and convert them into SMILES for data processing and analysis. 
A lot of these tools are well summarized by the Blue Obelisk Society Open Source Review [OBoyle:2016-9]. And they work to some degree of accuracy. 
These tools are then improved on and machine learning starts dominating as a model that sits on top to fix any inaccuracies of the algrothim. 
If we took it another direction, where we selectively aggregate data based on popularity, valuability over time, and organized to a degree of functionality but that much expertise amongst one person is not enough. You need many opinions to come to a standard set. 

The problem is with previous platforms was the lack of participants involved in any one project and a method to contribute back for it to grow. 
Most software and especially old software can be difficult to install and handle on top of modern technology thus driving the
need for something sustainable that can naturally grow. The chemical universe is large and too big for one person to fathom. 
It takes a multitude of chemical diversity expertise to put together a well-thought chemical list of most relative compounds to their respetive community.
To implement our idea we needed to pick a coding language that has the ability to write easy objects similar to a common language, English, for everyone to understand; Python.

<p align="center">
  <img width="1000" height="800" src="images/figures/figure_2.png">
  <i>Figure 1: Language Construction </i>
</p>

We also chose python because of it's distribution infrastructure to easily install objects installed on the cloud on PyPi. This acts as the "Printing Service" that enabled the mass distribution of the Gideon Bible and `GlobalChem` will follow in the same manner. 

# Methodology and Implementation

## Paper Selection Philosophy

Within academia, professors, post-doctorates, and graduate students, by nature of our work are required to read extensively about 
selective specific scientific fields. This in turn gives us an expert opinion in what data we value most. To start a thin layer data organization 
we begin by forming connections of most relevant data according to chemicals subfields. This is in accordance to the authorship
where each expertise opinion is recognized for different fields. A graph overview of the 
e layout in `GlobalChem`.

<p align="center">
  <img width="1000" height="1000" src="images/figures/figure_1.png">
  <i>Figure 2: Network Graph of GlobalChem</i>
</p>

## Object-Oriented Design

`GlobalChem` follows a simple object-oriented design where directories are the parent nodes and articles/books are leaf nodes. 
. In `Figure 2`, each leaf node is labeled appropriately as a class name to the paper it was referenced from. Each paper object
has either the functional groups that correspond to that paper's overall functionality in IUPAC, Preferred Name, Acronyms, SMILES, SMARTS
format. The choice for this design was that as more users contribute then can expand into different directories, add their own, 
and provide their respective popular chemical list. Each paper that is elected is converted into a `namespace` module an object
whose name is indicative of it's functionality. An example for the drug design community is the paper "Rings In Drugs" [Taylor:2014-6] whose
pythonic object equivalent is now "RingsInDrugs" with two functional methods that retreive the IUPAC:SMILES/SMARTS that was embedded
into the paper. 

## Manual SMILES abstraction

Papers are selected based on interested and relevance in the scientific community dictated by us, the authors, in respect to our field. 
The SMILES is abstracted in a variety of methods:

-  For simplistic molecules one representation of the SMILES can be directly translated using visual 
inspection. This worked for compounds usually at the beginning of a reported list that were small molecules.  

- For more complex molecules, the image can be redrawn into the free version of ChemDraw and then translated into a non-canonical version of SMILES. 

- For papers where the SMILES are written and the IUPAC is not known. We translate the SMILES into ChemDraw and then retrieve the name. 
Note that some of the names were modified based on human inspection in favor as well for preferred names. 

- For polymer papers, the site points were omitted from the name and some of the nomenclature adjusted for preferred names
over traditional. For example: 'yl' to mark site points for polymer connections was removed in favour of reduced english complexity. 

- Some SMILES were adjusted from their radical complement as they served as connection points. Some decisions were made to keep the radical component especially in the case if the IUPAC blue book common substituents. 
- 
- SMARTS strings were adapted from the SMILES using RDKit [@Landrum:2019-5]

# Data

At the time of writing the list now the list stands at:

| Chemical List                       | Languages                    | # of Entries | References               |
|-------------------------------------|------------------------------|--------------|--------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | 20           | Common Knowledge         |
| Essential Vitamins                  | Preferred Name/SMILES/SMARTS | 13           | Common Knowledge         |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | 42           | [Fulmer:2010-5]          |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | 94           | [OpenSmiles]             |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | 333          | [CRC:2004]               |
| Rings in Drugs                      | IUPAC/SMILES/SMARTS          | 92           | [Taylor:2014-6]          |
| Phase 2 Hetereocyclic Rings         | IUPAC/SMILES/SMARTS          | 19           | [Broughton:2004-9]       |
| Privileged Scaffolds                | IUPAC/SMILES/SMARTS          | 47           | [Welsch:2010-6]          |
| Common Warheads Covalent Inhibitors | IUPAC/SMILES/SMARTS          | 29           | [Gehringer:2019-6]       |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | 78           | [Hiorns:2019-6]          |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | 499          | [Takeuchi:2021-9]        |
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | 24           | [Petri:2020-12]          |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | 29           | [Hu:2021-3]              |
| BRaf Inhibitors                     | IUPAC/SMILES/SMARTS          | 54           | [Agianian:2018-6]        |
| Common Amino Acid Protecting Groups | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | [Isidro-Llobet:2009-6]   |
| Emerging Perfluoroalkyls            | IUPAC/SMILES/SMARTS          | 27           | [Pelch:2019-9]           |
| Chemicals For Clay Adsorption       | IUPAC/SMILES/SMARTS          | 33           | [Orr:2019-9]             |
| Common Regex Patterns               | Mol2                         | 1            |                          |

<p align="center">
  <i>Table 1: GlobalChem Object List</i>
</p>

# Tests & Applications

A total collection of 2153 IUPAC/Preferred Name/Acronym to SMILES/SMARTS was collected (with redundacy) and dispersed across 17 objects in
an organized fashion by subject. The code was refactored extensively to allow for ease of object addition according to subject
and functionality. To test the tolerance of these lists to other software we test on a couple of open source platforms to determine 
data interoperability. Although, it can be suggested that some of the software implemented should be expanded to perhaps
include functional groups that couldn't be parsed. 

## Cheminformatics Test

Two open-source cheminformatic platforms have taken staple as foundational tools: RDKit and Indigo. To test each SMILES string, each string gets
passed into a `Mol` RDKit object and `Indigo.loadMolecule()` object where any failures are recorded logged in Table 2. Cheminformatic interoperability
between different platforms promotes a unification. This can be expanded into OpenBabel and many others as a tolerance checker.

| Software | Number of Failed Compounds | Example Failed SMILES                                                                                                                                                                        |
|----------|----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| RDKit    | 11                         | 'CSi(C(C)(C)C)C', 'CC(Si(C1=CC=CC=C1)C2=CC=CC=C2)(C)C', 'n1cncn1', 'C&1&1&1',  'c1ccccc1C&1&1', 'n1cccc1', 'C&1&1&1&1', 'n1ccnc1C'                                                           |
| Indigo   | 8                          | 'C&1&1&1&1', 'C&1&1&1', 'c1ccccc1C&1&1', 'CSi(C(C)(C)C)C', 'CC(Si(C1=CC=CC=C1)C2=CC=CC=C2)(C)C', 'CSi(C(C)(C)C)C', 'CC(Si(C1=CC=CC=C1)C2=CC=CC=C2)(C)C'                                      |

<p align="center">
  <i>Table 2: GlobalChem Tolerable Results</i>
</p>


## ForceFields Test

Common chemical groups by nature should be of interest for development of force-fields for molecular dynamic simulations that
are trying to simulate chemicals and biological systems movements, velocity, and charge, which can also serve as dual interoperable test 
for SMILES strings. Popular force-fields such as General Amber ForceField (GAFF) [Wang:2004-7], Optimized Potentials for Liquid Simulations (OPLS)
[Jorgensen:1988-7], and Charmm General ForceField (CGenFF) [Vanommeslaeghe:2010-3] collect common chemical lists of use to 
construct the best chemical representative compounds of a particular space and estimate the geometry and charge of atoms
using their atom-typing engine. Due to accessibility, we focused primiarily on CGenFF to check it's tolerance level. We 
developed a private version of CGenFF that can estimate the atom types from `SDF` bond type column. This enabled us to pass the SMILES strings
through `RDKit` (by nature anything failed in RDKit fails in CGenFF) and transform to `SDF` to a `CGenFF` stream output. 


CGenFF was founded on drug-like molecules and to test it's capability of handling what's reported in literature it's performance
is tested in accordance with it's penalty score distribution. The lower the distribution is to 0 the more performant the 
forcefield is. The distributions are reported in accordance with bonds, angles, dihedrals, charge classifications of the charmm potential energy
equation. 

<p align="center">
  <img width="1400" height="400" src="images/figures/figure_3.png">
  <i>Figure 3: Charmm Potential Energy Geometry and Charge Classifications</i>
</p>

We looked at BRAF Kinases Inhibitors for Cancer (54), Privileged Scaffolds (47), Common Warheads (29), Emerging PerfluoroAlkyls (27).
Any kinase inhibitors should exhibit drug-like features similar to what was chosen to CGenFF, privileged scaffolds are any 
elected scaffolding produced by nature, warheads designed for covalent inhibition, and a stretch into herbicides and toxicity
that could be used to kill us. 3 Failures were captured primarily due to the sulphonyl `SH` and the salt `O` to potassium.

<p align="center">
  <img width="500" height="100" src="images/figures/figure_6.png">
  <br>
  <i>Figure 4: Failed CGenFF Compounds</i>
</p>

### Discussion 

From results suggested from `GlobalChem` we can suggest for cheminformatic toolkits to expand more on Silicon based datasets as
well as handle the ampersand `&` operator for materials. Diamond is a common carbon substance that is indicated on the OpenSMILES
documentation as a `C&1&1&1&1`, as the parsers suggest this seems to be a problem in both `RDKit` and `Indigo`. 

CGenFF is paired with a penalty scoring method for how wrong the software is predicting the geometry of a small molecule
[Vanommeslaeghe:2012-12], using the penalty score of different classifications we can suggest precision avenues for what to 
fix for the ForceField Development. By passing `GlobalChem` into `CGenFF` we can recommend compounds of usefulness for 
manual parameterization [Kumar:2020] without having to rely on a brute force approach. 

<p align="center">
  <img width="1000" height="650" src="images/figures/figure_5.png">
  <i>Figure 5: Penalty Score distributions</i>
</p>

We can suggest things like targeting one of the Emerging Perfluoroalkyls, owed to their primarily high penalty scores as 
a unique addition into the forcefield that is impactful to the environmental chemical hazard community. 

# Conclusion

`GlobalChem` serves a purpose of documenting what is common and relevant to different chemical communities in a distributable
easy format with objects classified as primary paper functionality with methods containing the chemical list that accodomates
said functionality. 

`GlobalChem` has several proposed applications into machine learning and artifical intelligence for drug pipelines as a classification layer,
it can help construct a cheminformatic analysis of functional groups on a chemical dataset without any knowledge being known prior, and lastly
it has a potential educational use in teaching functional groups and SMILES to any potential chemistry students. 

# Acknowledgements

Thank you to Jacob Weiner, Tyree Wilson, and Paul Shapiro for their helpful discussions into the usability and functionality of GlobalChem.
Thank you to Blake Printy, and Robert Zeigler for their influence for python objected-oriented design and distribution.
Thank you to the University of Maryland School of Pharmacy department for promoting a collaborative and useful space for 
academics. 
