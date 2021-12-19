---
title: 'Global-Chem: Collections of common small molecules and their SMILES/SMARTS to support diverse chemical communities'
tags:
  - Python
  - Cheminformatics
authors:
  - name: Suliman Sharif,  Elena Yi Chow, Asuka Orr, Aarion Romany, Aziza Frank, Shaoqi Zhan, Ruibin Liu, Sunhwan Jo, Alexander D. MacKerell Jr. 
    orcid: 0000-0002-1342-9258
    affiliation: 1
affiliations:
 - name: University of Maryland Baltimore, School of Pharmacy
   index: 1
date: 12/08/2021
bibliography: paper.bib
---

# Introduction

The in silico chemical universe is expanding rapidly as open access titan databases (Enamine Database (20 Billion) [@Gorgulla:2020-4],
Zinc Database (2 Billion) [Irwin:2020-12], PubMed Database (68 Million) [Roberts:2001-2]) and cheminformatic tools
to process, manipulate, and derive new compound structures are established. While this chemical data big bang has yielded ultra-large datasets they are based on ambiguous classification systems making it difficult to systematically organize them for specific uses 
[I don't understand what you are saying here] .

Previously, partial organizational attempts were made on PubMed, filling chemical data linkages for computational toxicology called Actor for a specific
refactored and refined effort [Judson:2019-9]. For the EnamineDB, a scaffold associated with biological activity was designed to target 
Toll-Like Receptors in an object-oriented fashion [Perez-Regidor:2016-9]. These organizational methods are difficult
to extend to other systems as well as can be difficult to implement given the large amount of data. In addition, the information content of these papers is of limited utility to the common developer. 

To organize data we apply the idea of communication. Humans use symbols and drawings to communicate, a set of symbols and the rules to combining them are called a language. 
are called a language. Different languages can be employed to carry different features and mean different things to a variety of communities as the infer meaning. 
IUPAC was a written language that predates even drawing atoms as a method of communication between chemists [Cooke-Fox:1989-5]; 
other chemical sub-communities also adopted the language and applied to their field to different dialects i.e polymer chemistry, organo-metallic chemistry.
In the recent years, SMILES [Weininger:1988-5] is becoming a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D or 3D geometry with ease.
Due to it's "first to market" scientific chemical language IUPAC is the legacy language that is a lexical key to unlocking informational wealth about a chemical pattern or group. But there are problems with the language due to it's length in describing bigger molecules. IUPAC names in organic chemisty papers can extend pages with no real value. To compact information, chemists just released the drawings but that can be hard to store precisely. Algorithms
are being designed to abstract and interpolate skeletal patterns and languages and convert them into SMILES for data processing and analysis. 
A lot of these tools are well summarized by the Blue Obelisk Society Open Source Review [OBoyle:2016-9]. And they work to some degree of accuracy. 
These tools are then improved on and machine learning starts dominating as a model that sits on top to fix any inaccuracies of the algrothim. 
If we took it another direction, where we selectively aggregate data based on popularity, valuability over time, and organized to a degree of functionality but that much expertise amongst one person is not enough.
You need many opinions to come to a standard set. 

In the context of a well-classified chemical database the major challenge is the enormity of the chemical universe. Accordingly, it takes a range of chemical expertise to put together a well-thought chemical list of compounds relevant to their respetive community. Thus, it is necessary for a large number of participants to contribute in order for such a database to grow. 
However, most software and especially old software can be difficult to install and handle on top of modern technology thus driving the
need for something sustainable that is readily accessible to potential participants, allowing the database to naturally grow.
This need motivated the development of the presented `Global-Chem` database

To implement `Global-Chem` we needed to pick a coding language that has the ability to write easy objects for particpants to understand; Python.

<p align="center">
  <img width="1000" height="800" src="images/figures/figure_2.png">
  <i>Figure 1: Language Construction </i>
</p>

Python was also chosen because of it's distribution infrastructure that allows for easy installation of objects available on the cloud. This 
allows `Global-Chem` to function as a free service behaving in the same manner of communication and mass distribution as the Gideon Bible. 

# Methodology and Implementation

## Paper Selection Philosophy

Scientists, by nature of thier work, are required to read extensively about 
selected scientific fields as well as access the associated data. This allows for scientists to develop expert knowledge in the fields and data they value most. This requires a thin layer data organization that allows for the relevant information and data to be readily accessed.
To achieve this we begin by forming connections of the most relevant data according to chemicals subfields that have been authored
by experts in the different fields. A graph overview of the Module layout in `Global-Chem`.

<p align="center">
  <img width="1000" height="1000" src="images/figures/figure_1.png">
  <i>Figure 2: Network Graph of Global-Chem</i>
</p>

## Object-Oriented Design

[Would it be possible to cross reference between datasets?  If you select one set of molecules, it would be nice to identify if any of those compounds are in other datasets; this could ultimately extended to cross references sub-strings]

`Global-Chem` follows a simple object-oriented design where directories are the parent nodes and articles or books are leaf nodes. In `Figure 2`, each leaf node is labeled appropriately as a class name to the reference-source paper or book. Each reference object
has either the functional groups that correspond to that paper's overall functionality in IUPAC, Preferred Name, Acronyms, SMILES, SMARTS
format. The motivation for this design was that as more users contribute they can expand into different directories, add their own directory, 
and provide their chemical list of interest. Each paper that is selected is converted into a `namespace` module, an object
whose name is indicative of it's functionality. An example for the drug design community is the paper "Rings In Drugs" [Taylor:2014-6] whose
python object equivalent is now "RingsInDrugs" with two functional methods that retrieve the  IUPAC:SMILES/SMARTS dictionary that was embedded included in the master object `Global-Chem`. 

## Manual SMILES abstraction

References and associatied compound lists are selected based on the interests of the contributing authors.  This should include consideration of relevance to the scientific community. 
The SMILES strings are abstracted in a variety of methods:

-  For simple molecules one representation of the SMILES can be directly translated using visual 
inspection. This is typically appropriate for compounds at the beginning of a reported list that were the most common denominator rings. 

- For complex molecules the image can be redrawn in the free version of ChemDraw and then translated into SMILES. 

- For sources where the SMILES are written and the IUPAC is not known the SMILES are translated into ChemDraw and the name retrieved. 
Note that some of the names may be modified based on human inspection in favor as well for preferred names. 

- For polymer papers, the site points were omitted from the name and some of the nomenclature adjusted for preferred names
over traditional. For example: 'yl' to mark site points for polymer connections was removed in favour of reduced english complexity. 

- Some SMILES were adjusted from their radical complement as they served as connection points. Some decisions were made to keep the radical component especially in the case if the IUPAC blue book common substituents. 
over traditional names. 

- SMARTS strings were adapted from the SMILES using RDKit [@Landrum:2019-5]

# Data

At the time of writing the list of objects includes:

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
| BRAF Inhibitors                     | IUPAC/SMILES/SMARTS          | 54           | [Agianian:2018-6]        |
| Common Amino Acid Protecting Groups | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | [Isidro-Llobet:2009-6]   |
| Emerging Perfluoroalkyls            | IUPAC/SMILES/SMARTS          | 27           | [Pelch:2019-9]           |
| Chemicals For Clay Adsorption       | IUPAC/SMILES/SMARTS          | 33           | [Orr:2019-9]             |
| Common Regex Patterns               | Mol2                         | 1            |                          |

<p align="center">
  <i>Table 1: GlobalChem Object List</i>
</p>

# Tests & Applications

A total collection of 2153 IUPAC/Preferred Name/Acronym to SMILES/SMARTS was collected (with redundacy) across 17 objects in
an organized fashion by subject. The code was refactored extensively to allow for ease of object addition according to subject and functionality.

## Results 

To test the utility of these lists with other software tests were performed on three open source platforms to determine 
data interoperability. Although, it can be suggested that some of the software implemented should be expanded to perhaps
include functional groups that couldn't be parsed. 

## Cheminformatics Test

Two open-source cheminformatic platforms our now widely considered as foundational tools: RDKit and Indigo. To test each SMILES string, each string gets
passed into a `Mol` RDKit object and `Indigo.loadMolecule()` object where any failures are recorded. Results on the number of failed compounds out of the 2153 compounds along with example of failed molecules is presented in Table 2. Cheminformatic interoperability
between different platforms promotes wider utilization. For example, OpenBabel is another utility that may used as a tolerance checker.

| Software | Number of Failed Compounds | Example of Failed SMILES                                |
|----------|----------------------------|---------------------------------------------------------|
| RDKit    | 11                         | 'CSi(C(C)(C)C)C', 'C&1&1&1&1',                          |
| Indigo   | 8                          | 'C&1&1&1&1', 'CC(Si(C1=CC=CC=C1)C2=CC=CC=C2)(C)C'       |

<p align="center">
  <i>Table 2: GlobalChem Tolerable Results</i>
</p>

## ForceFields Test

Access to broad collections of chemical groups will be of interest for development of force fields for molecular modeling and molecular dynamic simulations, allowing for studies on a wider range of chemicals and biological systems. The ability of a force field to treat molecules in the database can also serve as dual interoperable test 
for SMILES strings. Popular force fields such as General Amber ForceField (GAFF) [Wang:2004-7], Optimized Potentials for Liquid Simulations (OPLS)
[Jorgensen:1988-7], and Charmm General Force Field (CGenFF) [Vanommeslaeghe:2010-3] are based on collections of chemicals that are representative of the particular region of chemical space that the force field was designed to cover. 
In practice, this involves the atom-typing engine of each force field being applied to each molecule followed by assignment of the appropriate parameters.
This is to a large exten associated with the coverage of a force field.  Thus, the compound lists in Global-Chem can be used to identify specific regions of chemical space that have limited coverage and, therefore, represents future regions of chemical space for force field development. 
In the present study we used CGenFF to check it's tolerance level for the range of molecules currently in Global-Chem. To facilitate this an in-house extension of CGenFF was used that can assign atom types from `SDF` bond type column. This enabled us to pass the SMILES strings
through `RDKit` and transform `SDF` to a `CGenFF` stream output. 
The resulting failures are also presented in Table 2. It should be noted by nature of the data processing workflow the anything that fails in RDKit fails in CGenFF.
[need to add details for the individual lists, which could be pretty interesting; perhaps include this in the discussion]


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
that could be used to kill us. 3 Failures were captured primarily due to the sulphonyl `SH` and the activated oxygen to it's potassium salt counter anion in the perfluoroalkyl. 

<p align="center">
  <img width="500" height="100" src="images/figures/figure_6.png">
  <br>
  <i>Figure 4: Failed CGenFF Compounds</i>
</p>

### Discussion 

[You need to discuss the utility of Global-Chem and how that utility will increase as users add more chemical data.  And then present some example use cases, such as CGenFF]

From results suggested from `GlobalChem` we can suggest for cheminformatic toolkits to expand more on Silicon based datasets as
well as 

An interesting observation from the present analysis is the ability of tools to handle the ampersand `&` operator in SMILES for materials. For example, diamond is a common carbon substance whose SMILES strings is indicated in the OpenSMILES
documentation as a `C&1&1&1&1`. As shown in Table 2, this fails in both `RDKit` and `Indigo` indicating that improve handling of the `&` operator is required.

CGenFF includes a penalty scoring method based on the extent of analog new parameters have to those explicitly available in the force field [Dude, this is totally wrong: "how wrong the software is predicting the geometry of a small molecule"]
[Vanommeslaeghe:2012-12]. Accordingly, the penalty scores of different chemical classes can suggest those in need of force field development. Accordingly. passing `GlobalChem` into `CGenFF` can facilitate the identification of compounds for manual parameterization [Kumar:2020]. 

[If include a figure you need to discuss it!]

<p align="center">
  <img width="1000" height="750" src="images/figures/figure_5.png">
  <i>Figure 5: Penalty Score distributions</i>
</p>

We can suggest things like targeting one of the Emerging Perfluoroalkyls, owed to their primarily high penalty scores as 
a unique addition into the forcefield that is impactful to the environmental chemical hazard community. 

# Conclusion

`Global-Chem` serves the purpose of facilitating collecting, documenting and accessing different chemical communities as dictated by user input. It involves a distributable
easy format with objects classified as primary paper functionality with methods containing the chemical list that accodomates
said functionality. With respect to broader applicatibility `Global-Chem` will potentially be of utility for machine learning and artifical intelligence tools in drug development pipelines as a classification layer.  In addition, it can help construct a cheminformatic analysis of functional groups on a chemical dataset and, lastly,
it has a potential educational use in teaching functional groups and SMILES to any potential chemistry students. 

# Acknowledgements

Thank you to Jacob Weiner, Tyree Wilson, Paul Shapiro for their helpful discussions into the usability and functionality of Global-Chem.
Appreciation to the University of Maryland School of Pharmacy Depatment of Pharmaceutical Chemistry for promoting a collaborative and useful space for 
academics. Financial support from the NIH (GM131710) is acknowledged.

# Conflict of Interets

ADM is cofounder and CSO and SJ is Commercial Development Director of SilcsBio LLC.
