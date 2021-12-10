---
title: 'Global-Chem: A record collection of small molecules and their SMILES/SMARTS'
tags:
  - Python
  - Cheminformatics
authors:
  - name: Suliman Sharif, Shaoqi Zhan, Ruibin Liu, Sunhwan Jo, Aarion Romney, Elena Yi Chow, Alexander D. MacKerell Jr. 
    orcid: 0000-0002-1342-9258
    affiliation: 1
affiliations:
 - name: University of Maryland, School of Pharmacy
   index: 1
date: 12/08/2021
bibliography: paper.bib
---

# Introduction

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
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | kinase_electrophilic_warheads_smiles, kinase_electrophilic_warheads_smarts                                                 | 24           | [Petri:2020-12]          |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | kinase_privileged_scaffolds_smiles, kinase_privileged_scaffolds_smarts                                                     | 29           | [Hu:2021-3]              |
| BRaf Inhibitors                     | IUPAC/SMILES/SMARTS          | ribose_subpocket_smiles, ribose_subpocket_smarts, adenine_subpocket_smiles, adenine_subpocket_smarts                       | 20           | [Agianian:2018-6]        |
| BRaf Inhibitors (Continued)         | IUPAC/SMILES/SMARTS          | hydrophobic_subpocket_smiles, hydrophobic_subpocket_smarts,  type_1_subpocket_smiles, type_1_subpocket_smarts              | 17           |                          |
| BRaf Inhibitors (Continued)         | IUPAC/SMILES/SMARTS          | type_2_subpocket_smiles, type_2_subpocket_smarts, exposed_to_solvent_smiles, exposed_to_solvent_smarts                     | 17           |                          |
| Common Regex Patterns               | Mol2                         | common_regex_patterns                                                                                                      | 1            |                          |

# Conclusion

A total collection of 2534 key value pairs are recorded across 25 variables that can be accessed through the `GlobalChem` class. 
