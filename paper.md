---
title: 'Global-Chem: Collections of common small molecules and their SMILES/SMARTS to support diverse chemical communities'
tags:
  - Python
  - Cheminformatics
authors:
  - name: Suliman Sharif [first author]
    orcid: 0000-0002-1342-9258
    affiliation: 1
  - name: Elena Yi Chow
    orcid: 0000-0001-5559-1635
    affiliation: 2
  - name: Asuka Orr
    orcid: 0000-0003-4628-526X
    affiliation: 1
  - name: Aarion Romany
    affiliation: 1
  - name: Aziza Frank
    affiliation: 1
  - name: Shaoqi Zhan
    orcid: 0000-0002-6383-1771
    affiliation: 3
  - name: Ruibin Liu
    orcid: 0000-0001-8395-9353
    affiliation: 1
  - name: Sunhwan Jo
    affiliation: 2
  - name: Chris Burke
    affiliation: 2
  - name: Alexander D. MacKerell Jr. [corresponding author]
    affiliation: 1
affiliations:
  - name: University of Maryland Baltimore, School of Pharmacy
    index: 1
  - name: Independent Researcher
    index: 2
  - name: University of Oxford
    index: 3
date: 14 January 2022
bibliography: paper.bib
---

# Introduction

The in silico chemical universe is expanding rapidly as open access titan databases (Enamine Database (20 Billion) [@Gorgulla:2020-4],
Zinc Database (2 Billion) [@Irwin:2020-12], PubMed Database (68 Million) [@Roberts:2001-2]) and cheminformatic tools
to process, manipulate, and derive new compound structures are established. While this chemical data big bang has yielded useful ultra-large datasets they are based on ambiguous classification systems making it difficult to systematically organize them for specific uses.

<p align="center">
  <img width="400" height="300" src="images/figures/figure_4.png">
  <br />
  <i>Figure 1: Screenshot of the ZincDB request URLS</i>
</p>

For example, in `Figure 1`, the directory setup for downloading ZincDB molecules is shown. As is evident, the information content of the directory nomenclature does not contain information on the compounds they contain, making it nearly impossible to access specific molecules or classes molecules.  
Towards overcoming this, partial organizational attempts were made in PubMed, filling chemical data linkages for computational toxicology called Actor for a specific
refactored and refined effort [@Judson:2019-9]. In another example, for the EnamineDB a scaffold associated with biological activity was designed to target 
Toll-Like Receptors in an object-oriented fashion [@Perez-Regidor:2016-9]. However, these organizational methods are difficult
to extend to other systems and can be difficult to implement given the large amount of data.
In addition, the information content of these papers is of limited utility to the common developer. 

To organize chemical compounds we apply the idea of communication. Humans use symbols and drawings to communicate, a set of symbols and the rules to combining them are called a language. Languages can be employed to carry relevant, distinct features and mean something to their respective community. diagrammatically shown in `Figure 2`. 
International Union of Pure and Applied Chemistry (IUPAC) was a coalition that formed in the 1800s and their method of communication was named after the organization, IUPAC. 
IUPAC is a written language that predates even drawing atoms as a method of communication between chemists [@Cooke-Fox:1989-5]. 
Other chemical sub-communities adopted the IUPAC language and applied it to their fields that are comprised of different dialects i.e polymer chemistry, organo-metallic chemistry.
Due to it's "first to market" status, the scientific chemical language IUPAC is the legacy language that is the lexical key to unlocking information about a chemical pattern or group. 
But there are problems with the language due to it's length in describing bigger molecules. Simply, IUPAC names in organic chemisty papers are impractical, effecting extending the length of a manuscript,  while being of limited value given the challenge of interpreting such names.

To compact information, chemists presented drawings of chemical structures but information in such a format is hard to store precisely. Alternatively, SMILES [@Weininger:1988-5] has become a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D chemical connectivity information with ease.  Algorithms
have been designed to abstract and interpolate skeletal patterns and languages from chemical drawings and convert them into SMILES for data processing and analysis. 
A number of these tools, which work to varying degrees of accuracy, have been well summarized by the Blue Obelisk Society Open Source Review [@OBoyle:2016-9]. 
Efforts to improve these tools recently have included machine learning (ML) methods that essential "sit" on top of the underlying algoritm to fix any inaccuracies of the method. 
As an alternative we can take another direction, where data is selectively aggregated based on known classifications, popularity and utility, being organized to a degree of functionality that facilitates more widespread use. However, the criteria for such an aggregation of data is built upon human expertise, requiring input from a variety of people to attain the broadness and accessibility that would facilitate scientific discovery. 
In other words, in the context of a well-classified chemical database the major challenge is the enormity of the chemical universe, requring a range of chemical expertise to put togethe well-thought chemical lists of compounds relevant to their respective communities. Thus, it is necessary to create a tool to allow for a large number of participants to contribute in order for such a data compilation to grow. 
However, most software and especially old software can be difficult to install and handle on top of modern technology thus hindering participation. This situation drives the
need for a tool that is sustainable and readily accessible to potential participants, allowing the database to naturally grow.
This need motivated the development of the presented `GlobalChem` database tool.

To implement `GlobalChem` we selected a coding language that has the ability to write easy objects for particpants to understand; Python [10.5555/159351][@Cooke:1989-5].

<p align="center">
  <img width="700" height="450" src="images/figures/figure_2.png"><br/>
  <i>Figure 2: Language organized by category and functionality </i>
</p>
<br />
Python was also chosen because of it's distribution infrastructure that allows for easy installation of objects available on the cloud. This 
allows `Global-Chem` to function as a highly accessible tool that will allow users to readily access the chemical lists as well as to add content thereby continuosly expanding its utility. 

# Methodology and Implementation

## Chemical Set Selection & Object-Oriented Design Philosophy

Scientists, by nature of their work, are required to read extensively about 
selected scientific fields as well as access the associated data. This allows for scientists to develop expert knowledge in the fields and data they value most.
To take advantage of this knowledges requires a thin layer data organization that allows for the relevant information and data to be readily accessed.
To achieve this we begin by forming connections of the most relevant data according to chemicals sub-fields that have been authored
by experts in the different fields. `Figure 3` depicts the node Module layout of `Global-Chem`.  The layout shows an unweighted, 
arbitrary node hierarchy of the chemical sets included in `Global-Chem` as defined by the experts that introduce the data. Each blue circle represents a relevant field and their subsequent tree networks are highlighted by a contrasting colour.

<p align="center">
  <img width="1000" height="700" src="images/figures/figure_1_new.png">
  <i>Figure 3: Node Network of Global-Chem</i>
</p>

The tree network follows a simple object-oriented pythonic design in conjunction with literature where head nodes are the major corresponding scientific field (example: "Medicinal Chemistry") and their corresponding child nodes are the manuals, articles or books that are the references for the lists.
Each reference object has either the functional groups that correspond to that paper's overall functionality in IUPAC, Preferred Name, Acronyms, SMILES, or SMARTS
format. The motivation for this design was that as more users contribute they can expand into different directories, add their own directory, 
and provide their chemical list of interest. Each paper that is submitted is converted into a `namespace` module, an object
whose name is indicative of it's functionality. An example for the drug design community is the paper "Rings In Drugs" [@Taylor:2014-6] whose
python object equivalent is now "RingsInDrugs" with two functional methods that retrieve the  IUPAC:SMILES/SMARTS dictionary that was embedded included in the master object `Global-Chem`. 
Users can choose to cross reference leaf nodes between each other and do comparative chemical list studies since the IUPAC name and SMILES name are consistent across lists.
Note that not all the SMILES being portrayed are canonical given that users can create their own SMILES, which are not unique. To account for this users can parse `Global-Chem` SMILES into the `RDKit` parser
for canonical SMILES conversion. 

## Data

At the time of writing the list of nodes include those shown in Table 1. The list range from well defined classes of chemicals, such as amino acids, to more diverse lists such as Rings in Drugs. In addition, the languages used for each list are given, along with the number entires in the list and the reference. 
In addition, the number of times that compounds in each list fail in the CGenFF program, as discussed below, is given.

| Node List                           | Languages                    | # of Entries | References                |  CGenFF Errors            |
|-------------------------------------|------------------------------|--------------|---------------------------| --------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | 20           | Common Knowledge          | 0                         |
| Essential Vitamins                  | Preferred Name/SMILES/SMARTS | 13           | Common Knowledge          | 0                         |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | 42           | [@Fulmer:2010-5]          | 3                         |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | 94           | [@OpenSmiles]             | 10                        |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | 333          | [@CRC:2004]               | 1 (Excluding Radicals)    |
| Rings in Drugs                      | IUPAC/SMILES/SMARTS          | 92           | [@Taylor:2014-6]          | 0                         |
| Phase 2 Hetereocyclic Rings         | IUPAC/SMILES/SMARTS          | 19           | [@Broughton:2004-9]       | 0                         |
| Privileged Scaffolds                | IUPAC/SMILES/SMARTS          | 47           | [@Welsch:2010-6]          | 0                         |
| Common Warheads Covalent Inhibitors | IUPAC/SMILES/SMARTS          | 29           | [@Gehringer:2019-6]       | 4                         |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | 78           | [@Hiorns:2019-6]          | 7                         |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | 499          | [@Takeuchi:2021-9]        | 15                        |
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | 24           | [@Petri:2020-12]          | 0                         |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | 29           | [@Hu:2021-3]              | 0                         |
| BRAF Inhibitors                     | IUPAC/SMILES/SMARTS          | 54           | [@Agianian:2018-6]        | 5                         |
| Common Amino Acid Protecting Groups | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | [@Isidro-Llobet:2009-6]   | 41                        |
| Emerging Perfluoroalkyls            | IUPAC/SMILES/SMARTS          | 27           | [@Pelch:2019-9]           | 1                         |
| Chemicals For Clay Adsorption       | IUPAC/SMILES/SMARTS          | 33           | [@Orr:2019-9]             | 0                         |
| Schedule 1 United States Narcotics  | Preferred Name/SMILES/SMARTS | 240          | [@21CFRPart1]             | 1                         |
| Schedule 2 United States Narcotics  | Preferred Name/SMILES/SMARTS | 60           | [@21CFRPart1]             | 1                         |
| Schedule 3 United States Narcotics  | Preferred Name/SMILES/SMARTS | 22           | [@21CFRPart1]             | 1                         |
| Schedule 4 United States Narcotics  | Preferred Name/SMILES/SMARTS | 77           | [@21CFRPart1]             | 0                         |
| Schedule 5 United States Narcotics  | Preferred Name/SMILES/SMARTS | 8            | [@21CFRPart1]             | 0                         |
| Common Regex Patterns               | Mol2                         | 1            |                           | N/A                       |

<p align="center">
  <i>Table 1: GlobalChem Master Node Network</i>
</p>

## GlobalChemExtensions

To exhibit the wide functionality and uses cases of `GlobalChem` we created an extension cheminformatics tool independent of the Graph Network 
because it depends on other open source dependencies. The reasoning behind this was that sometimes users would just want
the data and to ease the installation process we concomitantly two different components that work together. Full exhibition of the
`GlobalChemExtensions` can be found in the `Gitbook` documentation available here (https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/). 
To highlight one functionality is the deep graph network extended into the plotly parallel coordinates plot shown in `Figure 8`. This gives
visibility into deep lexical layered graphs and help aid in organizing sets of chemical data.

<p align="center">
  <img width="1000" height="450" src="images/figures/figure_9.png">
  <br>
  <i>Figure 9: Plotly Conversion using `GlobalChemExtensions` </i>
</p>


The head node is `GlobalChem` and each subsequent layer is a "deep layer" that serves as nodes for a network. Users can build
their own networks and organize data as they see fit in a neural fashion. This helps expand chemical architectural neural 
strategies for node construction. Users can access the scatter deep layer functionality with [@Plotly] and others (Radial Analysis,
Principal Component Analysis, Language Conversion, Software Inteoperable Conversion between python objects) 

# Conclusion

`Global-Chem` was developed to facilitate accessing lists of known chemical compounds as objects to allow them to be used in the context of python-based workflows.
However, it can also facilitate the evaluation of other tools to access chemical information in the form of SMILES. An interesting observation from the present data is the ability of tools to handle the ampersand `&` operator in SMILES for materials. 
For example, diamond is a common carbon substance whose SMILES strings is indicated in the OpenSMILES
documentation as a `C&1&1&1&1`. As shown in Table 2, this fails in both `RDKit` and `Indigo` indicating that improved handling of the `&` operator is required. 

Beyond accessing SMILES stings we've shown the utility of `Global-Chem` to interogate the coverage of the force field `CGenFF`. By partitioning chemical space into well-defined chemical lists, `Global-Chem` allows for regions of chemical space where the CGenFF programs fails or assigns parameters of low analogy to be readily identified. This information will allow for decisions to be made concerning the addition of molecules in the CGenFF training set thereby allowing for systematic improvements in the force field.

## Statement of Need

`Global-Chem` was developed to facilitate the ability of scientists in both academia and industry to make their compounds of interest readily available to the scientific community in the form of objects that may be directly accessed from python. 
Accordingly, `Global-Chem` has a number of potential purposes, including teaching and cheminformatics, but our main perogative is to create a free record collection.
As `Global-Chem` requires direct user input, if we plant the seed now then, hopefully, our tree will grow. 
The actual growth of the tree will be decided on by the common chemical community and experts in the field. Enjoy. 

# Acknowledgements

Thank you to Daniel Khavrutskii, Jacob Weiner, Tyree Wilson, and Paul Shapiro for their helpful discussions into the usability and functionality of Global-Chem.
Appreciation to past mentors James Ryan, Robert Zeigler, and Blake Printy for discussions on good manufacturing practices of python packaging and distribution.
Appreciation to the University of Maryland Baltimore, School of Pharmacy, Department of Pharmaceutical Chemistry for promoting a collaborative and useful space for 
academics. Financial support from the NIH (GM131710) is acknowledged.

# Conflict of Interets

ADM is cofounder and CSO and SJ is Commercial Development Director of SilcsBio LLC. Chris Burke is Senior DevOps Engineer at L7 Informatics. 
