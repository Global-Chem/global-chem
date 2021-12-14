---
title: 'Global-Chem: A record collection of common small molecules and their SMILES/SMARTS in different chemical communities'
tags:
  - Python
  - Cheminformatics
authors:
  - name: Suliman Sharif,  Elena Yi Chow, Sunhwan Jo, Shaoqi Zhan, Ruibin Liu, Aarion Romany, Aziza Frank, Asuka Orr, Alexander D. MacKerell Jr. 
    orcid: 0000-0002-1342-9258
    affiliation: 1
affiliations:
 - name: University of Maryland, School of Pharmacy
   index: 1
date: 12/08/2021
bibliography: paper.bib
---

# Introduction

The chemical universe is expanding rapidly as open access titan databases (Enamine Database (20 Billion) [@Gorgulla:2020-4],
Zinc Database (2 Billion) [Irwin:2020-12], PubMed Database (68 Million) [Roberts:2001-2]) and cheminformatic tools
to process, manipulate, and derive new compound structures are established. This left us with a chemical data big bang
with ultra-large datasets and an ambiguous classification system in an attempt to organize the data. Previously, partial
organizational attempts were made on PubMed filling chemical data linkages for computational toxicology called Actor for a specific
refactored and refined effort [Judson:2019-9]. For the EnamineDB, a scaffold to biological activity was designed to target 
Toll-Like Receptors in an object-oriented fashion [Perez-Regidor:2016-9]. These organizational methods are difficult
to reproduce as well as can be difficult to implement given the amount of data. When applying these papers they don't provide
so much use to the common developer. So what do we do?

To organize the data we need to revert back to the idea of communication. Humans use symbols and drawings to communicate, a collection of symbols and their combinations
are called a language. Different languages can be employed to carry different features and mean different things to a variety of communities. 
IUPAC was a written language that predates even drawing atoms as a method of communication between chemists [Cooke-Fox:1989-5]; 
other chemical sub-communities also adopted the language and applied to their field to different dialects i.e polymer chemistry, organo-metallic chemistry.
In the recent years, SMILES [Weininger:1988-5] is becoming a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D or 3D geometry with ease.
Unfortunately, IUPAC is a legacy language and is the lexical key to informational wealth about a chemical pattern or group. Algorithms
were designed to abstract and interpolate skeletal patterns and languages and convert them into SMILES for data processing and analysis. 
A lot of these tools are well summarized by the Blue Obelisk Society Open Source Review [OBoyle:2016-9].

The problem is the lack of participants involved in any one project and a method to contribute back for it to grow. 
Most software and especially old software can be difficult to install and handle on top of modern technology thus driving the
need for something sustainable that can naturally grow. The chemical universe is large and too big for one person to fathom. 
It takes a multitude of chemical diversity expertise to put together a well-thought chemical list of most relative compounds to their respetive community.
To implement our idea we needed to pick a coding language that has the ability to write easy objects for everyone to understand; Python.

<p align="center">
  <img width="1000" height="800" src="images/figures/figure_2.png">
  <i>Figure 1: Language Construction </i>
</p>

We also chose python because of it's distribution infrastructure to easily install objects installed on the cloud. This 
acts a free service where `GlobalChem` will behave in the same manner as the Gideon Bible. 

# Methodology and Implementation

## Paper Selection Philosophy

Within academia, professors, post-doctorates, and graduate students, by nature of our work are required to read extensively about 
selective specific scientific fields. This in turn gives us an expert opinion in what data we value most. To start a thin layer data organization 
we begin by forming connections of most relevant data according to chemicals subfields. This is in accordance to the authorship
where each expertise opinion is recognized for different fields. A graph overview of the Module layout in `GlobalChem`.

<p align="center">
  <img width="800" height="800" src="images/figures/figure_1.png">
  <i>Figure 2: Network Graph of GlobalChem</i>
</p>

## Object-Oriented Design

`GlobalChem` follows a simple object-oriented design where directories are the parent nodes and articles/books are leaf nodes. 
. In `Figure 2`, each leaf node is labeled appropiately as a class name to the paper it was referenced from. Each paper object
has either the functional groups that correspond to that paper's overall functionality in IUPAC, Preferred Name, Acronyms, SMILES, SMARTS
format. The choicefor this design was that as more users contribute then can expand into different directories, add their own, and provide their own chemical list.

### Data abstraction

Papers are selected based on interested and relevance in the scientific community dictated by us, the authors, in respect to our field. 
The SMILES is abstracted in a variety of methods:

-  For simplistic molecules one representation of the SMILES can be directly translated using visual 
inspection. This worked for compounds usually at the beginning of a reported list that were the most common denominator rings. 

- For complex molecules, the image can be redrawn into the free version of ChemDraw and then translated into SMILES. 

- For papers where the SMILES are written and the IUPAC is not known. We translate the SMILES into ChemDraw and then retrieve the name. 
Note that some of the names were modified based on human inspection in favor as well for preferred names. 

- For polymer papers, the site points were omitted from the name and some of the nomenclature adjusted for preferred names
over traditional. 

- SMARTS strings were adapted from the SMILES using RDKit [@Landrum:2019-5]

# Results

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
| Common Warheads                     | IUPAC/SMILES/SMARTS          | 29           | [Gehringer:2019-6]       |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | 78           | [Hiorns:2019-6]          |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | 499          | [Takeuchi:2021-9]        |
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | 24           | [Petri:2020-12]          |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | 29           | [Hu:2021-3]              |
| BRaf Inhibitors                     | IUPAC/SMILES/SMARTS          | 54           | [Agianian:2018-6]        |
| Common Amino Acid Protecting Groups | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | [Isidro-Llobet:2009-6]   |
| Emerging Perfluoroalkyls            | IUPAC/SMILES/SMARTS          | 27           | [Pelch:2019-9]           |
| Chemicals For Clay Adsorption       | IUPAC/SMILES/SMARTS          | 33           | [Orr:2019-9]             |
| Common Regex Patterns               | Mol2                         | 1            |                          |

# Conclusion

A total collection of 2534 key value pairs are recorded across 25 variables that can be accessed through the `GlobalChem` class. 
