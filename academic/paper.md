---
title: 'Global-Chem: A General Public Dictionary of Common Chemical Names to SMILES'

authors:
  - name: Suliman Sharif [first author]
    orcid: 0000-0002-1342-9258
    affiliation: 1
  - name: Ruibin Liu
    orcid: 0000-0001-8395-9353
    affiliation: 1
  - name: Asuka A. Orr
    orcid: 0000-0003-4628-526X
    affiliation: 1
  - name: Daniel Khavrutskii
    orcid: 0000-0001-5014-5572
    affiliation: 7
  - name: Sunhwan Jo
    affiliation: 2
  - name: Bettina Lier
    orcid: 0000-0002-8032-0084
    affiliation: 4
  - name: Chris Burke
    affiliation: 2
  - name: Jacob Weiner
    orcid: 0000-0003-0033-6033
    affiliation: 1
  - name: Anastasia Croitoru
    orcid: XXX
    affiliation: 5
  - name: Aziza Frank
    affiliation: 1
  - name: Nathaniel McClean
    orcid: XXX
    affiliation: 1
  - name: Anmol Kumar
    orcid: XXX
    affiliation: 1
  - name: Aarion Romany
    affiliation: 1
  - name: Mingtian Zhao
    orcid: XXX
    affiliation: 1
  - name: Takayuki Serizawa
    affiliation: 2
  - name: Jared Deacon
    affiliation: 8
  - name: Ian Jones
    affiliation: 1
  - name: Shaoqi Zhan
    orcid: 0000-0002-6383-1771
    affiliation: 3
  - name: Mike Woster
    affiliation: 6
  - name: Rebecca Pinette-Dorin
    affiliation: 6
  - name: Elena Yi Chow
    orcid: 0000-0001-5559-1635
    affiliation: 2
  - name: Sevien Schulhoff
    orcid: XXXX
    affiliation: 1
  - name: Alexander D. MacKerell Jr. [corresponding author]
    orcid: 0000-0001-8287-6804
    affiliation: 1
affiliations:
  - name: University of Maryland Baltimore, School of Pharmacy
    index: 1
  - name: Independent Researcher
    index: 2
  - name: University of Oxford
    index: 3
  - name: University of Natural Resources and Life Sciences, Vienna
    index: 4
  - name: École Polytechnique
    index: 5
  - name: Linux Foundation
    index: 6
  - name: University of Maryland, Baltimore County
    index: 7
  - name: University of Maryland, College Park
    index: 8
---

# Introduction

The *in silico* chemical universe is expanding rapidly as open access titan databases such as Enamine Database (20 Billion) (1),
Zinc Database (2 Billion) (2), and PubMed Database (68 Million) (3), as well as cheminformatics tools
are being developed for processing, manipulating, and deriving new compound structures. This big bang of chemical data has produced useful ultra-large data sets, but based on ambiguous classification systems that make it difficult to organize them systematically and select molecules of interest based on the name. Previously developed organizational methods are hard to distribute to a mass audience and complicated to implement given a large amount of data. In addition, the information content is of limited use because the software infrastructure on which it was developed is difficult for the ordinary developer to reproduce. 

A resource for the selection of chemicals that are of "general" purpose is needed by various sub-communities to relate data acros fields. International Union of Pure and Applied Chemistry (IUPAC) is a written language that predates even drawing atoms as a method of communication between chemists (7). Over time, the IUPAC names naturally turned into a slang ("preferred") language due to humans wanting to speak it while communicating with each other. Effectively, the **natural** chemical language that is extant today is a blend of both formal and informal nomenclature with no defined rules. To visualize chemical information, chemists presented drawings of chemical structures (skeletal diagrams). However, in recent modern times, drawings are not easy to query on a mass scale. Alternatively, Simplified Molecular Input Line Entry Specification (SMILES, 8) was initiated for the organic chemist as a different way to encode compounds on the computer. SMILES strings were designed with a set of grammar rules and became a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D chemical connectivity information. Selecting chemical compounds for a given purpose requires expertise. Expertise is gained by experience and studying a dedicated discipline. Dedicated disciplines most often have a set of common functional groups that are relevant to that community. This allows us to focus on valuable compounds. We do not need all the compounds, since a lot of them are not useful or not possible. In our paper, we describe how `Global-Chem`, an open source chemical dictionary , was developed to facilitate the ability of scientists in both academia and industry to make their compounds of interest readily available to the scientific community in the form of objects that may be directly accessed from python. 

<p align="center">
<img width="1000" alt="Screen Shot 2022-07-16 at 5 29 41 PM" src="https://user-images.githubusercontent.com/11812946/179372511-61758864-6b0a-410e-b15f-578fd8227a14.png">
    <br>
    <i>Figure 1: File Architecture of Global-Chem</i>
</p>

# Data & Features

At the time of writing the list of objects includes those shown in Table 1. The list ranges from well defined natural classes of chemicals, such as amino acids, vitamins or salts to more diverse lists such as rings in drugs, emerging perfluoroalkyls etc. In addition, the languages used for each list are given, along with the number of entries in the list and the reference. The number of times that compounds in each list fail in our automated SMILES to CGenFF workflow.

| Chemical List                        | Languages                    | # of Entries | References               |  
|--------------------------------------|------------------------------|--------------|--------------------------| 
| Amino Acids                          | IUPAC/SMILES/SMARTS          | 20           | Common Knowledge         | 
| Essential Vitamins                   | Preferred Name/SMILES/SMARTS | 13           | Common Knowledge         | 
| Common Organic Solvents              | IUPAC/SMILES/SMARTS          | 42           | (13)                     | 
| OPEN SMILES                          | IUPAC/SMILES/SMARTS          | 94           | (14)                     | 
| IUPAC Blue Book (CRC Handbook) 2003  | Preferred Name/SMILES/SMARTS | 333          | (15)                     | 
| Rings in Drugs                       | IUPAC/SMILES/SMARTS          | 92           | (16)                     | 
| Phase 2 Heterocyclic Rings           | IUPAC/SMILES/SMARTS          | 19           | (17)                     | 
| Privileged Scaffolds                 | IUPAC/SMILES/SMARTS          | 47           | (18)                     | 
| Common Warheads Covalent Inhibitors  | IUPAC/SMILES/SMARTS          | 29           | (19)                     | 
| Common Polymer Repeating Units       | IUPAC/SMILES/SMARTS          | 78           | (20)                     | 
| Common R Group Replacements          | IUPAC/SMILES/SMARTS          | 499          | (21)                     | 
| Electrophilic Warheads for Kinases   | Preferred Name/SMILES/SMARTS | 24           | (22)                     | 
| Privileged Scaffolds for Kinases     | IUPAC/SMILES/SMARTS          | 29           | (23)                     | 
| BRAF Inhibitors                      | IUPAC/SMILES/SMARTS          | 54           | (24)                     | 
| Common Amino Acid Protecting Groups  | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | (25)                     | 
| Emerging Perfluoroalkyls             | IUPAC/SMILES/SMARTS          | 27           | (26)                     | 
| Chemicals For Clay Adsorption        | IUPAC/SMILES/SMARTS          | 33           | (27)                     | 
| Schedule 1 United States Narcotics   | Preferred Name/SMILES/SMARTS | 240          | (28)                     | 
| Schedule 2 United States Narcotics   | Preferred Name/SMILES/SMARTS | 60           | (28)                     | 
| Schedule 3 United States Narcotics   | Preferred Name/SMILES/SMARTS | 22           | (28)                     | 
| Schedule 4 United States Narcotics   | Preferred Name/SMILES/SMARTS | 77           | (28)                     | 
| Schedule 5 United States Narcotics   | Preferred Name/SMILES/SMARTS | 8            | (28)                     | 
| PihKal                               | Preferred Name/SMILES/SMARTS | 179          | (29)                     | 
| Excipients Cimetidine & Acyclovir    | Preferred Name/SMILES/SMARTS | 14           | (30)                     | 
| HowToLiveLonger	                     | Preferred Name/SMILES/SMARTS | 4            | (31)                     | 
| Monoclonal Antibodies                | Preferred Name/SMILES/SMARTS | 19           | (32)                     | 
| Common Lubricants for Sex Wellness   | Preferred Name/SMILES/SMARTS | 38           | (33)                     |
| FDA Tainted Sexual Enhancements      | Preferred Name/SMILES/SMARTS | 4            | (34)                     | 
| Common Food Salts                    | Preferred Name/SMILES/SMARTS | 14           | (35)                     | 
| FDA Color Additive List 1            | FDA Name/SMILES/SMARTS       | 12           | (36)                     | 
| FDA Color Additive List 2            | FDA Name/SMILES/SMARTS       | 15           | (36)                     |
| FDA Color Additive List 3            | FDA Name/SMILES/SMARTS       | 16           | (36)                     | 
| FDA Color Additive List 4            | FDA Name/SMILES/SMARTS       | 39           | (36)                     | 
| FDA Color Additive List 5            | FDA Name/SMILES/SMARTS       | 27           | (36)                     | 
| FDA Color Additive List 6            | FDA Name/SMILES/SMARTS       | 29           | (36)                     | 
| FDA Color Additive List 7            | FDA Name/SMILES/SMARTS       | 37           | (36)                     | 
| Constituents of Cannabis Sativa      | Name/SMILES/SMARTS           | 394          | (37)                     |
| Phytocannabinoids                    | Name/SMILES/SMARTS           | 111          | (38)                     |
| Organophosphorous Nerve Toxic Agents | Name/SMILES/SMARTS           | 14           | (39)                     |
| Cengage Bronsted Acids               | Name/SMILES/SMARTS           | 42           | (40)                     |
| Common Regex Patterns                | Mol2                         | 1            |                          |

<p align="center">
  <i>Table 1: Global-Chem Object List Columns: "Chemical List" is the name of the node that contains the chemical list, "Languages" specifies the name and their respective translations, "Number of Entries" is how many molecules exist within one node, "References" are the resources the molecules were recorded from, and the last column "CGenFF Errors" is how many times CGenFF skipped a molecule. If the value is "N/A" it means it was a node added after testing and allows room for additional chemical space exploration.</i>
</p>

At the time of writing the list of features includes those shown in Table 2. The list range from well defined algorithms implemented into Global-Chem and their respective description and discipline.

| Software Feature              | Description                                                                    |  Code Length | Discipline | Reference |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|---------- |-----------|--------------|
| Validating SMILES             | An adapter to other SMILES platforms (RDKit, PySMILES, SELFIES, PartialSMILES, DeepSMILES, MolVS) to validate by interoperability           | 107        | Cheminformatics          |  (41), (42), (43), (44), (45), (46) | 
| Decoding Fingerprints and Classifying SMILES            | Decoding fingerprints to complex SMILES and to an IUPAC using an annotated dictionary of bit vectors | 129        | Cheminformatics          | (47) |
| SMILES Bidirectional PDF Parsing      | Converting lists of SMILES to 2D drawings in PDF parsing and parsing PDF back to SMILES | 685        | Cheminformatics          | (48) |
| Drug Design Filtering      | Filtering lists of SMILES by a variety of common drug filters (Lipinski Rule of 5, Ghose, Veber, Rule of 3, REOS, Drug-Like, Filters | 137        | Cheminformatics           | (49), (50), (51), (52), (53), (54) |
| Deep Layer Scattering     | Scattering Nodes of Collective SMILES and their relations to each other in Parallel Coordinate Diagram implemented in Plotly | 184        | Cheminformatics          | (55) |
| SMARTS Identifier     | A web application implemented in Flask to test the validation of SMARTS submatching of strings against the MiniFrag Database | 307        | Cheminformatics          |  (56) |
| Protonating SMILES    | A distributable version of the Dimorphite-DL package to protonate SMILES over a range of pH with a control over variant production | 56        | Cheminformatics          | (57) |
| Sunbursting SMILES    | Applying a sunburst plot to large collection of SMILES to identify functional groups and pairs of functional groups within the set | 253        | Cheminformatics          | (58), (59) |
| Peptide Sequence to SMILES    | An evolution of Cocktail-Shaker to include Lanthipeptides and covalent sulphur linkages in SMILES strings | 147        | Cheminformatics     | (60) |
| Visualization SMARTS    | A python application programming interface to port the SMARTS.plus visualizer for SMARTS strings into a jupyter notebook | 47        | Cheminformatics          | (61) |
| One-Hot Encoding SMILES    | A Global-Chem encoder that encodes SMILES for Machine Learning including the '&' denoted as a polymer ex. Diamond | 112        | Cheminformatics          | (62) | 
| Principal Component Analysis on SMILES    | A principal component analysis on a list of SMILES with hyperparamter tuning for morgan fingerprinting provided and visualization with Bokeh | 154        | Cheminformatics          | (63), (64) |
| Networkx Adapter    | A graph to graph network adapter between Global-Chem and NetworkX for ease of interoperability for data engineering | 65        | Cheminformatics          | (65) |
| Scaffold Graph Adapter   | An adapter to take a large collection of Global-Chem Nodes and analyze their Structure Hierachy with Scaffold Graphs | 97        | Cheminformatics          | (66) |
| Global-Chem Protein  | An adapter to biopandas to process pdb protein files as well as an implementation of the Bostrom Algorithm to Structurally Filter SMILES | 467        | Bioinformatics          | (67), (68) |
| Global-Chem RNA  | Conversion of RNA Sequence to SMILES and a visualizer for RNA sequences for Python Jupyter Notebooks | 181        | Bioinformatics          | (69) |
| Global-Chem DNA  | Conversion of DNA Sequence to SMILES and a visualizer for DNA sequences for Python Jupyter Notebooks | 181        | Bioinformatics          | (69) | 
| Global-Chem Bacteria  | A python model with attributes for general bacteria classifications as well as a common list | 214        | Bioinformatics          | (70) |
| Global-Chem Monoclonal Antibodies  | A python model with attributes for general monoclonoal antibodies classifications as well as a common list | 20        | Bioinformatics          | (71) | 
| Z-Matrix Store  | A python model store where users can pull standard z-matrices for molecules queried by their IUPAC | 159        | Quantum Chemistry          | (72) |
| Psi4 Parser  | A python model for analyzing psi4 output files and plotting interaction energy data automatically through Plotly | 193        | Quantum Chemistry          | (72) |
| Moly Adapter  | A software adapter and enhanced functionality for Moly and visualizing HOMO/LUMO orbitals of molecules | 87        | Quantum Chemistry          | (73) |
| Global-Chem Molecule  | A Global Molecule that can parse SMILES, GAFF, CGenFF Stream files into Pandas dataframes, a visualizer with Atom Types and SMILES in RDKit, new mix of cross-discipling languages (SMILES and CGenFF Atom Types) using CXSMILES, CXSMARTS, and CurlySMILES | 386        | Force Fields         | (74), (75), (76) | 
| CGenFF Molecule  | A CGenFF Parser that can parse, write edit, and update stream files with Pandas DataFrames | 532        | Force Fields         | (77) |
| GAFF2 Molecule   | A GAFF2 Parser that can parse, write edit, and update stream files with Pandas DataFrames  | 454        | Force Fields         | (77) |
| CGenFF Disimiliarity Score   | A CGenFF dissimilarity algorithm based on the atom types and their tuples of bonded parameters (bonds, angles, dihedrals, impropers) to determine a dissimilarity score | 191        | Force Fields         | (78) |
| Open Source Database Monitor   | An open source database monitor that performs heartbeat checks on common chemical lists running on cloud web servers | 95        | Development Operations         | (79) |
| Plotly Templates   | A Graphing template to use for Plotly to make your data look "pretty" | 80        | Graphing         |

<p align="center">
  <i>Table 2: Global-Chem-Extensions Feature List Columns: "Feature" name of the feature model, "Description" a summarized account of what the feature does, "Feature Code Length" is how many lines the actual feature occupies without including infrastructure, "Discipline" is what scientific discipline and distribution pathway does the feature exist, and the last column "References" is what scientific resource, if any, does the feature stem from.</i>
</p>

# Performance

Global-Chem SMILES strings are only valid as they are interoperable with other open source cheminformatic software. In Table 3, we evaluate the validity of SMILES strings in different sets of standards of the cheminformatic field.

| Software                        | # of Errors |  Examples             |
|---------------------------------| ----------- | --------------------- |
| RDKit                           |             |                       |
| PySMILES                        |             |                       |
| SELFIES                         |             |                       |
| PartialSMILES                   |             |                       |
| DeepSMILES                      |             |                       |
| MolVS                           |             |                       |

# Conclusion

Accordingly, `Global-Chem` has several potential purposes but our main prerogative is to create a free record collection. As `Global-Chem` requires direct user input, if we plant the seed now then, hopefully, our tree will grow. The actual growth of the tree will be decided on by the experts of the community dedicated to their field. Please Enjoy a Printable Data Format of the Dictionary here:  .

# Acknowledgements

Thank you to Tyree Wilson, Garrick Centola and Paul Shapiro for their helpful discussions on the usability, functionality, and guidance of the manuscript of Global-Chem. Thank you to Ryland Forsythe for our discussions on polymers and for extending SMILES to materials in Global-Chem. Thank you to Payal Chatterjee for discussions of the 1,3-dithiolane parametrization. Thank you to Steven Fletcher for our discussion on aziridine during an organic synthesis lecture. Thank you to Holden Smith for helping in the usability of sunbursting chemical information. Thank you to Melissa Landolf, Raina Landolf, and Ella Landolf for discussions about the theory and application to useful products. Appreciation to past mentors James Ryan, Robert Zeigler, and Blake Printy for discussions on good manufacturing practices of python packaging and distribution. Appreciation to the University of Maryland Baltimore, School of Pharmacy, Department of Pharmaceutical Chemistry for promoting a collaborative and useful space for academics from all different scientific disciplines. Financial support from the NIH (GM131710) is acknowledged.

# References

(1) Gorgulla, C.; et. al. An Open-Source Drug Discovery Platform Enables Ultra-Large Virtual Screens. Nature 2020, 580 (7805), 663–668. https://doi.org/10.1038/s41586-020-2117-z.

(2) Irwin, John. J.; et. al. ZINC20—A Free Ultralarge-Scale Chemical Database for Ligand Discovery. ACS Publications 2013, 60 (12), 6065–6073. https://doi.org/10.1021/acs.jcim.0c00675.

(3) Roberts, R. J. PubMed Central: The GenBank of the Published Literature. Proceedings of the National Academy of Sciences of the United States of America 2001, 98 (2), 381–382. https://doi.org/10.1073/pnas.98.2.381.

(4) Wang, Junmei, et al. “Development and Testing of a General Amber Force Field.” Journal of Computational Chemistry, vol. 25, no. 9, July 2004, pp. 1157–74. PubMed

(5) Jorgensen, William L., and Julian Tirado-Rives. “The OPLS [Optimized Potentials for Liquid Simulations] Potential Functions for Proteins, Energy Minimizations for Crystals of Cyclic Peptides and Crambin.” Journal of the American Chemical Society, vol. 110, no. 6, Mar. 1988, pp. 1657–66.

(6) Vanommeslaeghe, K., E. et al. CHARMM General Force Field (CGenFF): A Force Field for Drug-like Molecules Compatible with the CHARMM All-Atom Additive Biological Force Fields.” Journal of Computational Chemistry, vol. 31, no. 4, Mar. 2010, pp. 671–90. PubMed Central

(7) Cooke-Fox, D. I.; et al. Computer Translation of IUPAC Systematic Organic Chemical Nomenclature. 1. Introduction and Background to a Grammar-Based Approach. Journal of Chemical Information and Computer Sciences 1989, 29 (2), 101–105. https://doi.org/10.1021/ci00062a009.

(8) Weininger, D.; et al. SMILES, a Chemical Language and Information System. 1. Introduction to Methodology and Encoding Rules. Journal of Chemical Information and Computer Sciences 1988, 28 (1), 31–36. https://doi.org/10.1021/ci00062a009.

(9) Vanommeslaeghe, K., and A. D. MacKerell. “Automation of the CHARMM General Force Field (CGenFF) I: Bond Perception and Atom Typing.” Journal of Chemical Information and Modeling, vol. 52, no. 12, Dec. 2012, pp. 3144–54

(10) Vanommeslaeghe, K., E. Prabhu Raman, et al. “Automation of the CHARMM General Force Field (CGenFF) II: Assignment of Bonded Parameters and Partial Atomic Charges.” Journal of Chemical Information and Modeling, vol. 52, no. 12, Dec. 2012, pp. 3155–68. 

(11) Sharif Suliman,  Jo Sunhwan, MacKerell Jr. Alexander D. "Automation of SMILES to CGenFF" Unpublished, 2021. 

(12) Mackerell, Alexander D. “Empirical Force Fields for Biological Macromolecules: Overview and Issues.” Journal of Computational Chemistry, vol. 25, no. 13, Oct. 2004, pp. 1584–604. PubMed,

(13) Fulmer, G. R.; et al. NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist. Organometallics 2010, 29 (9), 2176–2179. https://doi.org/10.1021/om100106e.

(14) Daylight Theory. OpenSmiles.

(15) Lide, D. R.; et al. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data; CRC Press, 2004. https://doi.org/10.5860/choice.37-4225.

(16) Taylor, R. D.; et al. Rings in Drugs. Journal of Medicinal Chemistry 2014, 57 (14), 5845–5859. https://doi.org/10.1021/jm4017625.

(17) Broughton, H. B.; Watson, I. A. Selection of Heterocycles for Drug Design.; PubMed, 2004; Vol. 23, pp 51–58. https://doi.org/10.1016/j.jmgm.2004.03.016.

(18) Welsch, M. E.; et al. Privileged Scaffolds for Library Design and Drug Discovery.; PubMed, 2010; Vol. 14, pp 347–361. https://doi.org/10.1016/j.cbpa.2010.02.018.

(19) Gehringer, M. E.; Laufer, S. A. Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology; ACS Publications, 2019; Vol. 62, pp 5673–5724. https://doi.org/10.1021/acs.jmedchem.8b01153.

(20) Hiorns, R. C. A Brief Guide to Polymer Nomenclature (IUPAC Technical Report).; 2012; Vol. 84, pp 2167–2169. https://doi.org/10.1351/PAC-REP-12-03-05.

(21) Takeuchi, K.; et al. R-Group Replacement Database for Medicinal Chemistry.; 2021; Vol. 7, p FSO742. https://doi.org/10.2144/fsoa-2021-0062.

(22) Petri, L.; et al. An Electrophilic Warhead Library for Mapping the Reactivity and Accessibility of Tractable Cysteines in Protein Kinases.; 2020; Vol. 207, p 112836. https://doi.org/10.1016/j.ejmech.2020.112836.

(23) Hu, H.; et al. Systematic Comparison of Competitive and Allosteric Kinase Inhibitors Reveals Common Structural Characteristics.; 2021; Vol. 214, p 113206. https://doi.org/10.1016/j.ejmech.2021.113206.

(24) Agianian, B.; Gavathiotis, E. Current Insights of BRAF Inhibitors in Cancer.; ACS Publications, 2018; Vol. 61, pp 5775–5793. https://doi.org/10.1021/acs.jmedchem.7b01306.

(25) Isidro-Llobet, A.; et al. Amino Acid-Protecting Groups.; ACS Publications, 2009; Vol. 109, pp 2455–2504. https://doi.org/10.1021/cr800323s.

(26) Pelch, K.; et al. PFAS Health Effects Database: Protocol for a Systematic Evidence Map.; Science Direct, 2019; Vol. 130, p 104851. https://doi.org/10.1016/j.envint.2019.05.045.

(27) Orr, A.; et al. Combining Experimental Isotherms, Minimalistic Simulations, and a Model to Understand and Predict Chemical Adsorption onto Montmorillonite Clays.; PubMed, 2021; Vol. 6, pp 14090–14103. https://doi.org/10.1021/acsomega.1c00481.

(28) ECFR :: 21 CFR Part 1308 - Schedules.

(29) Shulgin, Alexander T., and Ann Shulgin. Pihkal: A Chemical Love Story. 1. ed., 8. print, Transform, 2010.

(30) Vaithianathan, Soundarya, et al. “Effect of Common Excipients on the Oral Drug Absorption of Biopharmaceutics Classification System Class 3 Drugs Cimetidine and Acyclovir.” Journal of Pharmaceutical Sciences, vol. 105, no. 2, Feb. 2016, pp. 996–1005. PubMed

(31) "How to Live Longer" Github, https://github.com/geekan/HowToLiveLonger

(32) Data Abstracted from https://labels.fda.gov/

(33) EXSENS-USA.com. “Lube Lessons 3: The Sex Lube Ingredient Glossary.” EXSENS-USA.Com,

(34) Research, Center for Drug Evaluation and. “Tainted Sexual Enhancement Products.” FDA, June 2022

(35) Belot, Laure  "Alimentation : face aux doutes, les internautes s'organisent". Le Monde.

(36) https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list

(37) Turner, C. E., et al. “Constituents of Cannabis Sativa L. XVII. A Review of the Natural Constituents.” Journal of Natural Products, vol. 43, no. 2, Apr. 1980, pp. 169–234. PubMed

(38) Hanuš, Lumír Ondřej, et al. “Phytocannabinoids: A Unified Critical Inventory.” Natural Product Reports, vol. 33, no. 12, Nov. 2016, pp. 1357–92. PubMed,

(39) Mukherjee, Sudisha, and Rinkoo Devi Gupta. “Organophosphorus Nerve Agents: Types, Toxicity, and Treatments.” Journal of Toxicology, vol. 2020, Sept. 2020, p. 3007984.

(40) PKa Values for Organic and Inorganic Bronsted Acids at 25 Celsius.

(41) The RDKit 2021.03.1 Documentation

(42) "Pysmiles: The lightweight and pure-python SMILES reader and writer" Github, https://github.com/pckroon/pysmiles

(43) Krenn, Mario, et al. “Self-Referencing Embedded Strings (SELFIES): A 100% Robust Molecular String Representation.” Machine Learning: Science and Technology, vol. 1, no. 4, Dec. 2020, p. 045024.

(44) "Partialsmiles" Github, https://github.com/baoilleach/partialsmiles

(45) O’Boyle, Noel, and Andrew Dalke. DeepSMILES: An Adaptation of SMILES for Use in Machine-Learning of Chemical Structures. Sept. 2018. chemrxiv.org,

(46) "MolVS: Molecule Validation and Standardization." Github, https://github.com/mcs07/MolVS

(47) Rogers, David, and Mathew Hahn. “Extended-Connectivity Fingerprints.” Journal of Chemical Information and Modeling, vol. 50, no. 5, May 2010, pp. 742–54.

(48) "Epam/Indigo: Universal Cheminformatics Libraries, Utilities and Database Search Tools.” GitHub, https://github.com/epam/Indigo

(49) Lipinski, C. A., et al. “Experimental and Computational Approaches to Estimate Solubility and Permeability in Drug Discovery and Development Settings.” Advanced Drug Delivery Reviews, vol. 46, no. 1–3, Mar. 2001, pp. 3–26. PubMed

(50) Ghose, A. K., et al. “A Knowledge-Based Approach in Designing Combinatorial or Medicinal Chemistry Libraries for Drug Discovery. 1. A Qualitative and Quantitative Characterization of Known Drug Databases.” Journal of Combinatorial Chemistry, vol. 1, no. 1, Jan. 1999, pp. 55–68. PubMed

(51) Veber, Daniel F., et al. “Molecular Properties That Influence the Oral Bioavailability of Drug Candidates.” Journal of Medicinal Chemistry, vol. 45, no. 12, June 2002, pp. 2615–23. PubMed

(52) Congreve, Miles, et al. “A ‘Rule of Three’ for Fragment-Based Lead Discovery?” Drug Discovery Today, vol. 8, no. 19, Oct. 2003, pp. 876–77. 

(53) Walters, W. Patrick, and Mark Namchuk. “Designing Screens: How to Make Your Hits a Hit.” Nature Reviews Drug Discovery, vol. 2, no. 4, Apr. 2003, pp. 259–66.

(54) Bickerton, G. Richard, et al. “Quantifying the Chemical Beauty of Drugs.” Nature Chemistry, vol. 4, no. 2, Jan. 2012, pp. 90–98. PubMed Central.

(55) Wilkinson, Leland, and Graham Wills. The Grammar of Graphics. 2nd ed, Springer, 2005.

(56) O’Reilly, Marc, et al. “Crystallographic Screening Using Ultra-Low-Molecular-Weight Ligands to Guide Drug Design.” Drug Discovery Today, vol. 24, no. 5, May 2019, pp. 1081–86. ScienceDirect

(57) Ropp, Patrick J., et al. “Dimorphite-DL: An Open-Source Program for Enumerating the Ionization States of Drug-like Small Molecules.” Journal of Cheminformatics, vol. 11, no. 1, Feb. 2019, p. 14. BioMed Central

(58) Corp, (c) 2008-2018 Software Ambience. “DaisyDisk, the Most Popular Disk Analyzer for Mac.” DaisyDisk

(59) Sunburst. https://plotly.com/python/sunburst-charts/

(60) Sharif, Suliman. “Cocktail Shaker: An Open Source Drug Expansion and Enumeration Library for Peptides.” Journal of Open Source Software, vol. 5, no. 52, Aug. 2020, p. 1992. joss.theoj.org

(61) Ehrt, Christiane, et al. “SMARTS.plus – A Toolbox for Chemical Pattern Design.” Molecular Informatics, vol. 39, no. 12, Dec. 2020, p. 2000216

(62) Franky. “Basic Molecular Representation for Machine Learning.” Medium, 20 Sept. 2021,

(63) Ding, Chris, and Xiaofeng He. “K -Means Clustering via Principal Component Analysis.” Twenty-First International Conference on Machine Learning  - ICML ’04, ACM Press, 2004, p. 29. 

(64) A Step-by-Step Explanation of Principal Component Analysis (PCA) | Built In.

(65) Hagberg Aric A., Schult Daniel A., and Swart Pieter J., “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008

(66) Scott Oliver B., A W Edith Chan, "ScaffoldGraph: an open-source library for the generation and analysis of molecular scaffold networks and scaffold trees, Bioinformatics" Vol. 36, 12, 15 June 2020, pp. 3930–3931, https://doi.org/10.1093/bioinformatics/btaa219

(67) Raschka, Sebastian. “BioPandas: Working with Molecular Structures in Pandas DataFrames.” Journal of Open Source Software, vol. 2, no. 14, June 2017, p. 279. joss.theoj.org

(68) Boström, Jonas, et al. “Do Structurally Similar Ligands Bind in a Similar Fashion?” Journal of Medicinal Chemistry, vol. 49, no. 23, Nov. 2006, pp. 6716–25.

(69) "DNA Features Viewer" Github, https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer

(70) Pitt, T. L., and M. R. Barer. “Classification, Identification and Typing of Micro-Organisms.” Medical Microbiology, 2012, pp. 24–38. PubMed Central

(71) Orr, Asuka, Sharif Suliman "Common Monoclonal Antibodies", Unpublished. 

(72) Turney, Justin M., et al. “Psi4: An Open-Source Ab Initio Electronic Structure Program: Psi4: An Electronic Structure Program.” Wiley Interdisciplinary Reviews: Computational Molecular Science, vol. 2, no. 4, July 2012, pp. 556–65. 

(73) "Moly - Molecular Visualization in Jupyter" Github, https://github.com/VHchavez/moly

(74) Chemaxon Extended SMILES and SMARTS - CXSMILES and CXSMARTS | Chemaxon Docs

(75) Drefahl, Axel. “CurlySMILES: A Chemical Language to Customize and Annotate Encodings of Molecular and Nanodevice Structures.” Journal of Cheminformatics, vol. 3, Jan. 2011, p. 1. PubMed Central

(76) Barr, Avron, and Edward A. Feigenbaum. The Handbook of Artificial Intelligence. Volume I Volume I. 1981.

(77) Reback, Jeff, et al. Pandas-Dev/Pandas: Pandas 1.4.3. v1.4.3, Zenodo, 2022.

(78) Hall, Lowell H., et al. “Molecular Similarity Based on Novel Atom-Type Electrotopological State Indices.” Journal of Chemical Information and Computer Sciences, vol. 35, no. 6, Nov. 1995, pp. 1074–80

(79) "Upptime" Github, https://github.com/upptime/upptime

(80) Chemical Weapons Disarmament in Russia: Problems and Prospects. Henry L. Stimson Center, 1995.

(81) Kumar, Anmol, et al. “FFParam: Standalone Package for CHARMM Additive and Drude Polarizable Force Field Parametrization of Small Molecules.” Journal of Computational Chemistry, vol. 41, no. 9, 2020, pp. 958–70. Wiley Online Library

(82) "About FFParam-GUI", Ffparam 1.0 Documentation.

(83) Gaussian 16, Revision C.01, Frisch, M. J.; Trucks, G. W.; Schlegel, H. B.; Scuseria, G. E.; Robb, M. A.; Cheeseman, J. R.; Scalmani, G.; Barone, V.; Petersson, G. A.; Nakatsuji, H.; Li, X.; Caricato, M.; Marenich, A. V.; Bloino, J.; Janesko, B. G.; Gomperts, R.; Mennucci, B.; Hratchian, H. P.; Ortiz, J. V.; Izmaylov, A. F.; Sonnenberg, J. L.; Williams-Young, D.; Ding, F.; Lipparini, F.; Egidi, F.; Goings, J.; Peng, B.; Petrone, A.; Henderson, T.; Ranasinghe, D.; Zakrzewski, V. G.; Gao, J.; Rega, N.; Zheng, G.; Liang, W.; Hada, M.; Ehara, M.; Toyota, K.; Fukuda, R.; Hasegawa, J.; Ishida, M.; Nakajima, T.; Honda, Y.; Kitao, O.; Nakai, H.; Vreven, T.; Throssell, K.; Montgomery, J. A., Jr.; Peralta, J. E.; Ogliaro, F.; Bearpark, M. J.; Heyd, J. J.; Brothers, E. N.; Kudin, K. N.; Staroverov, V. N.; Keith, T. A.; Kobayashi, R.; Normand, J.; Raghavachari, K.; Rendell, A. P.; Burant, J. C.; Iyengar, S. S.; Tomasi, J.; Cossi, M.; Millam, J. M.; Klene, M.; Adamo, C.; Cammi, R.; Ochterski, J. W.; Martin, R. L.; Morokuma, K.; Farkas, O.; Foresman, J. B.; Fox, D. J. Gaussian, Inc., Wallingford CT, 2016.

(84) Møller, Chr., and M. S. Plesset. “Note on an Approximation Treatment for Many-Electron Systems.” Physical Review, vol. 46, no. 7, Oct. 1934, pp. 618–22. 

(85) Dunning, Thom H. “Gaussian Basis Sets for Use in Correlated Molecular Calculations. I. The Atoms Boron through Neon and Hydrogen.” The Journal of Chemical Physics, vol. 90, no. 2, Jan. 1989, pp. 1007–23

(86) Brooks, B. R., et al. “CHARMM: The Biomolecular Simulation Program.” Journal of Computational Chemistry, vol. 30, no. 10, July 2009, pp. 1545–614. PubMed Central

(87) Vanommeslaeghe, Kenno, et al. “Robustness in the Fitting of Molecular Mechanics Parameters.” Journal of Computational Chemistry, vol. 36, no. 14, May 2015, pp. 1083–101. PubMed

(88) Guvench, Olgun, and Alexander D. MacKerell. “Automated Conformational Energy Fitting for Force-Field Development.” Journal of Molecular Modeling, vol. 14, no. 8, Aug. 2008, pp. 667–79. Springer Link

(89) Draghici, Cristian, and Jon T. Njardarson. “Chemistry By Design: A Web-Based Educational Flashcard for Exploring Synthetic Organic Chemistry.” Journal of Chemical Education, vol. 89, no. 8, July 2012, pp. 1080–82.

(90) Morgan, H. L. “The Generation of a Unique Machine Description for Chemical Structures-A Technique Developed at Chemical Abstracts Service.” Journal of Chemical Documentation, vol. 5, no. 2, May 1965, pp. 107–13

(91) Levenshtein, V. I. “Binary Codes Capable of Correcting Deletions, Insertions and Reversals.” Soviet Physics Doklady, vol. 10, Feb. 1966, p. 707. 

(92) Mobley, David L., et al. “Escaping Atom Types in Force Fields Using Direct Chemical Perception.” Journal of Chemical Theory and Computation, vol. 14, no. 11, Nov. 2018, pp. 6076–92.

(93) Tol, Paul. 2021. “Colour Schemes.” Technical note SRON/EPS/TN/09-002 3.2. SRON.

(94) Mozilla Public License, Version 2.0. https://www.mozilla.org/en-US/MPL/2.0/.

(95) Commissioner, Office of the. “Part 11, Electronic Records; Electronic Signatures - Scope and Application.” U.S. Food and Drug Administration, 11 June 2020,

(96) Rhodes, Colin, et al. “Regulatory Compliance Requirements for an Open Source Electronic Image Trial Management System.” Conference Proceedings : ... Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE Engineering in Medicine and Biology Society. Annual Conference, vol. 2010, 2010, pp. 3475–78. PubMed Centra

(97). Jolliffe, I. T., editor. “Principal Component Analysis and Factor Analysis.” Principal Component Analysis, Springer, 2002, pp. 150–66. Springer Link, https://doi.org/10.1007/0-387-22440-8_7.

(98) Ding, Chris, and Xiaofeng He. “K -Means Clustering via Principal Component Analysis.” Twenty-First International Conference on Machine Learning  - ICML ’04, ACM Press, 2004, p. 29.

(99) Leach, Paul J., et al. A Universally Unique IDentifier (UUID) URN Namespace. Request for Comments, RFC 4122, Internet Engineering Task Force, July 2005. IETF,

# Conflict of Interets

ADM is cofounder and CSO, and SJ is Commercial Development Director of SilcsBio LLC. Chris Burke is Senior DevOps Engineer at L7 Informatics. Mike Woster is the Chief Revenue Officer of the Linux Foundation, Rebecca Pinette-Dorin is marketing researcher at Exsens.
