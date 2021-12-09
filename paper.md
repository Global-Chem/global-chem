---
title: 'Global-Chem: A record collection of small molecules and their SMILES/SMARTS'
tags:
  - Python
  - Cheminformatics
authors:
  - name: Suliman Sharif, Shaoqi Zhan, Ruibin Liu, Sunhwan Jo, Elena Yi Chow, Alexander D. MacKerell Jr. 
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
Zinc Database (2 Billion) [Irwin:2020-12], PubMed Database (68 Million) [Roberts:2001-2]) and cheminformatic tools to process, manipulate, and derive new
compound structures are established. With these new methods we now have a combinatorial explosion leaving us with ultra-large datasets
and no proper organization. This presents a problem, with clustering and statistical methods being applied to chemical data and no actual way to organize
compounds into their respective functional groups and their functionality leaves little to no understanding and hard to make use of. 

To classify we need to revert back to the idea of communication. Humans use symbols and drawings to communicate, a collection of symbols and their combinations
are called a language. Different languages can be employede to carry different features and mean different things to a variety of communities. 
For organic chemistry, we draw skeletal patterns to communicate, and a written language using english characters to communicate, this was established as IUPAC [Cooke-Fox:1989-5],
other chemical subcommunities also adopted the language and applied to their field in different dialects i.e polymer chemistry, organometallic chemistry.
In the recent years, SMILES [Weininger:1988-5] is becoming a popular 1-D language amongst cheminformaticians for large datasets compared to 2D or 3D data (more lines). 
With most ultra-large datasets stored primarily in SMILES, we can start at organizing the data into appropiate classifiers that
will make it useful for all communities.

To accomplish this, we make use of on the most useful papers that classify a functionality of a set of compounds. If we 
write the IUPAC name and the SMILES name we can create a 1:1 mapping as a piece in the puzzle of the drug pipeline. The SMILES can be used
to connect to other resources in the chemical data space to physics. Whereas the IUPAC name will talk to the chemists and 
possibly some biologists. We then move it into a easy distributable format where it's accessible to everyone to expand.

# Methodology and Implementation

### Software

``` GlobalChem``` operates as a standalone python3 object and packaged distributed on PyPi. It has a no dependencies 
as it's just a key, value store of strings. Different papers are organized as retrieval functions on the object which fetch
classifier nomenclature (similar to the paper it was abstracted from) in IUPAC/SMILES and IUPAC/SMARTS. The SMARTS string
is a superset of SMILES and is used as a string matcher when looking for functional groups. It can be derived using a 
popular chem(o)informatic tool, RDKit [@Landrum:2019-5].

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

# Results

At the time of writing the list now the list stands at:

| Chemical List                       | Languages                    | Variables                                                                                                                  | # of Entries | References               |
|-------------------------------------|------------------------------|----------------------------------------------------------------------------------------------------------------------------|--------------|--------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | amino_acid_side_smiles, amino_acid_side_smarts                                                                             | 20           | Common Knowledge         |
| Essential Vitamins                  | Preferred Name/SMILES/SMARTS | vitamins_smiles, vitamins_smarts                                                                                           | 13           | Common Knowledge         |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | common_organic_solvents_smiles, common_organic_solvents_smarts                                                             | 42           | [Fulmer:2010-5]          |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | open_smiles_functional_groups_smiles, open_smiles_functional_groups_smarts                                                 | 94           | [OpenSmiles]             |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | iupac_blue_book_radical_smiles, iupac_blue_book_radical_smarts, iupac_blue_book_rings_smiles, iupac_blue_book_rings_smarts | 333          | [CRC:2004]               |
| Rings in Drugs                      | IUPAC/SMILES/SMARTS          | rings_in_drugs_smiles, rings_in_drugs_smarts                                                                               | 92           | [Taylor:2014-6]          |
| Phase 2 Hetereocyclic Rings         | IUPAC/SMILES/SMARTS          | phase_2_hetereocyclic_rings_smiles, phase_2_hetereocyclic_rings_smarts                                                     | 19           | [Broughton:2004-9]       |
| Privileged Scaffolds                | IUPAC/SMILES/SMARTS          | privileged_scaffolds_smiles, privileged_scaffolds_smarts                                                                   | 47           | [Welsch:2010-6]          |
| Common Warheads                     | IUPAC/SMILES/SMARTS          | common_warhead_smiles, common_warhead_smarts                                                                               | 29           | [Gehringer:2019-6]       |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | common_polymer_repeating_units_smiles, common_polymer_repeating_units_smarts                                               | 78           | [Hiorns:2019-6]          |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | r_groups_replacements_smiles, r_groups_replacements_smarts                                                                 | 499          | [Takeuchi:2021-9]        |
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | kinase_electrophilic_warheads_smiles, kinase_electrophilic_warheads_smarts                                                 | 29           | [Petri:2020-12]          |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | kinase_privileged_scaffolds_smiles, kinase_privileged_scaffolds_smarts                                                     | 29           | [Hu:2021-3]              |
| Common Regex Patterns               | Mol2                         | common_regex_patterns                                                                                                      | 1            |                          |

# Conclusion

A total collection of 2534 key value pairs are recorded across 25 variables that can be accessed through the `GlobalChem` class. 
