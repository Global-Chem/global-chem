---
title: 'Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities'

authors:
  - name: Suliman Sharif [first author]
    orcid: 0000-0002-1342-9258
    affiliation: 1
  - name: Ruibin Liu
    orcid: 0000-0001-8395-9353
    affiliation: 1
  - name: Asuka Autumn Grace Orr
    orcid: 0000-0003-4628-526X
    affiliation: 1
  - name: Daniel Khavrutskii
    orcid: XXX
    affiliation: 5
  - name: Bettina Lier
    affiliation: 4
  - name: Anastasia Croitoru
    orcid: XXX
    affiliation: 5
  - name: Chris Burke
    affiliation: 2
  - name: Aziza Frank
    affiliation: 1
  - name: Sunhwan Jo
    affiliation: 2
  - name: Jacob Weiner
    orcid: 0000-0003-0033-6033
    affiliation: 1
  - name: Nathaniel McClean
    orcid: XXX
    affiliation: 1
  - name: Aarion Romany
    affiliation: 1
  - name: Mingtian Zhao
    orcid: XXX
    affiliation: 1
  - name: Takayuki Serizawa
    affiliation: 6
  - name: Jared Deacon
    affiliation: 6
  - name: Ian Jones
    affiliation: 1
  - name: Shaoqi Zhan
    orcid: 0000-0002-6383-1771
    affiliation: 3
  - name: Anmol Kumar
    affiliation: 6
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
    affiliation: 1
affiliations:
  - name: University of Maryland Baltimore, School of Pharmacy
    index: 1
  - name: Independent Researcher
    index: 2
  - name: University of Oxford
    index: 3
  - name: University of Vienna
    index: 4
  - name: École Polytechnique
    index: 5
  - name: Takeda Pharmaceuticals
    index: 6
---

# Introduction

The *in silico* chemical universe is expanding rapidly as open access titan databases Enamine Database (20 Billion) (1),
Zinc Database (2 Billion) (2), PubMed Database (68 Million) (3) and cheminformatic tools
to process, manipulate, and derive new compound structures are established. While this chemical data big bang has yielded useful ultra-large datasets they are based on ambiguous classification systems making it difficult to systematically organize them and select molecules of interest based on the name. Previous organizational methods are hard to distribute to a mass audience and can be complicated to implement given the large amount of data. In addition, the information content of these papers is of limited utility because the infrastructure it was designed is hard to reproduce to the common developer. 

General Force Fields are designed to be "general purpose" force fields that capture chemical environments in the form of a dictionary of atom types and their physical properties in respect to themselves and each other. Major Force Fields such as General Amber Force Field (35), Optimized Potentials for Liquid Simulations (36), and CHARMM General Force Field (CGenFF, 37) were originally contrived from chemical lists that were "drug-like" since it was designed for the medicinal chemistry community. If we would like to expand CGenFF into moving beyond drugs we need to elect a variety of chemical space then molecules that serve as the best representation of their field. 

To select chemical compounds we apply the idea of communication. International Union of Pure and Applied Chemistry (IUPAC) is a written language that predates even drawing atoms as a method of communication between chemists (4). Over time, the IUPAC names naturally turned into a slang ("preferred") language due to humans wanting to speak it while communicating with each other. Effectively, the **natural** chemical language that is extant today is a blend of both a formal and informal nomenclature. To visualize chemical information, chemists presented drawings of chemical structures (skeletal diagrams) but information but, in recent modern times, drawings are not easy to query on a mass scale. Alternatively, a different way to encode compounds into the computer SMILES (5) designed for the organic chemist with a set of grammar rules became a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D chemical connectivity information. 

Penalty scores are attributed by the CGenFF program to molecules whose entire chemical connectivity is not present in CGenFF. When a molecule is passed to the CGenFF program it navigates through a set of rules that represent an atom type similarity network, a tree structure, similar to the Global-Chem file structure. Once its atom types along with chemical connectivity are known, bonded parameters available in CGenFF are assigned to the molecule. If exact matches of the bonded parameters are not available, a second tree traversal browses for alternate parameter by using a second-rules file that assigns penalties based on the analogy to known parameters. Once the bonded parameter substitutions with the lowest penalty score  are determined, the CGenFF Program assigns those parameters along with the associated penalties. In addition, the program identifies the original parameters that are also output into the stream file used in the various molecular modeling programs. Partial atomic charges and associated penalties are assigned through an extended bond-charge increment scheme. It consists in associating atom type along with chemical connectivity (including  bond, angle and dihedral) with charge increment values subtracted from the atoms formal charge. Thus, while the CGenFF program can successfully ingest a large number of molecules, the majority of those molecules are assigned penalties that indicate the level of analogy of the assigned bonded parameters and charges. 
 
To capture enough chemical lists we communicated, debated and elected papers of relevance and wrote a dictionary of IUPAC/Preferred to SMILES **manually** to effectively learn the language, organize the chemicals appropiately by their purpose, and develop rules to ensure data integrity and consistency. We developed an in-house automated workflow to process IUPAC/Preferred to SMILES to CGenFF (40) to enable us to safely exchange chemical information between collaborators by something we can naturally speak to select the most relevant compounds using the Penalty Score as a positive or negative feedback for forcefield paramitirization. 

Selecting chemical compounds requires expertise. Expertise is gained by experience and studying a dedicated discipline. Dedicated displines most often have a set of common functional groups that are relevant to that community, this allows us to focus on compounds that are valuable. We do not need all the compounds, since a lot of them are not useful or not possible. In our paper, we describe how `Global-Chem`, an open source knowledge graph, was developed to facilitate the ability of scientists in both academia and industry to make their compounds of interest readily available to the scientific community in the form of objects that may be directly accessed from python. However, research was directed under the guidance of the first author whose opinion may lead to some bias in the present methodolgy. 

<p align="center">
<img width="1000" alt="Screen Shot 2022-07-16 at 5 29 41 PM" src="https://user-images.githubusercontent.com/11812946/179372511-61758864-6b0a-410e-b15f-578fd8227a14.png">
    <br>
    <i>Figure 1: Knowledge Graph of GlobalChem</i>
</p>

# Data & Software Features

At the time of writing the list of objects include those shown in Table 1. The list range from well defined natural classes of chemicals, such as amino acids, vitamins, salt, and to more diverse lists such as rings in drugs, emerging perfluoroalkyls etc. In addition, the languages used for each list are given, along with the number entires in the list and the reference. The number of times that compounds in each list fail in our automated SMILES to CGenFF workflow.


| Chemical List                        | Languages                    | # of Entries | References               |  CGenFF Errors            |
|--------------------------------------|------------------------------|--------------|--------------------------| --------------------------|
| Amino Acids                          | IUPAC/SMILES/SMARTS          | 20           | Common Knowledge         | 0                         |
| Essential Vitamins                   | Preferred Name/SMILES/SMARTS | 13           | Common Knowledge         | 0                         |
| Common Organic Solvents              | IUPAC/SMILES/SMARTS          | 42           | (6)                      | 3                         |
| Open Smiles                          | IUPAC/SMILES/SMARTS          | 94           | (7)                      | 10                        |
| IUPAC Blue Book (CRC Handbook) 2003  | Preferred Name/SMILES/SMARTS | 333          | (8)                      | 1 (Excluding Radicals)    |
| Rings in Drugs                       | IUPAC/SMILES/SMARTS          | 92           | (9)                      | 0                         |
| Phase 2 Hetereocyclic Rings          | IUPAC/SMILES/SMARTS          | 19           | (10)                     | 0                         |
| Privileged Scaffolds                 | IUPAC/SMILES/SMARTS          | 47           | (11)                     | 0                         |
| Common Warheads Covalent Inhibitors  | IUPAC/SMILES/SMARTS          | 29           | (12)                     | 4                         |
| Common Polymer Repeating Units       | IUPAC/SMILES/SMARTS          | 78           | (13)                     | 7                         |
| Common R Group Replacements          | IUPAC/SMILES/SMARTS          | 499          | (14)                     | 15                        |
| Electrophillic Warheads for Kinases  | Preferred Name/SMILES/SMARTS | 24           | (15)                     | 0                         |
| Privileged Scaffolds for Kinases     | IUPAC/SMILES/SMARTS          | 29           | (16)                     | 0                         |
| BRAF Inhibitors                      | IUPAC/SMILES/SMARTS          | 54           | (17)                     | 5                         |
| Common Amino Acid Protecting Groups  | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | (18)                     | 41                        |
| Emerging Perfluoroalkyls             | IUPAC/SMILES/SMARTS          | 27           | (19)                     | 1                         |
| Chemicals For Clay Adsorption        | IUPAC/SMILES/SMARTS          | 33           | (20)                     | 0                         |
| Schedule 1 United States Narcotics   | Preferred Name/SMILES/SMARTS | 240          | (21)                     | 1                         |
| Schedule 2 United States Narcotics   | Preferred Name/SMILES/SMARTS | 60           | (21)                     | 1                         |
| Schedule 3 United States Narcotics   | Preferred Name/SMILES/SMARTS | 22           | (21)                     | 1                         |
| Schedule 4 United States Narcotics   | Preferred Name/SMILES/SMARTS | 77           | (21)                     | 0                         |
| Schedule 5 United States Narcotics   | Preferred Name/SMILES/SMARTS | 8            | (21)                     | 0                         |
| PihKal                               | Preferred Name/SMILES/SMARTS | 179          | (22)                     | 0                         |
| Excipients Cimetidine & Acyclovir    | Preferred Name/SMILES/SMARTS | 14           | (23)                     | N/A                       |
| HowToLiveLonger	                   | Preferred Name/SMILES/SMARTS | 4            | (24)                     | N/A                       |
| Monoclonal Antibodies                | Preferred Name/SMILES/SMARTS | 19           | (25)                     | N/A                       |
| Common Lubricants for Sex Wellness   | Preferred Name/SMILES/SMARTS | 38           | (26)                     | N/A                       |
| FDA Tainted Sexual Enhancements      | Preferred Name/SMILES/SMARTS | 4            | (27)                     | N/A                       |
| Common Food Salts                    | Preferred Name/SMILES/SMARTS | 14           | (28)                     | N/A                       |
| FDA Color Additive List 1            | FDA Name/SMILES/SMARTS       | 12           | (29)                     | N/A                       |
| FDA Color Additive List 2            | FDA Name/SMILES/SMARTS       | 15           | (29)                     | N/A                       |
| FDA Color Additive List 3            | FDA Name/SMILES/SMARTS       | 16           | (29)                     | N/A                       |
| FDA Color Additive List 4            | FDA Name/SMILES/SMARTS       | 39           | (29)                     | N/A                       |
| FDA Color Additive List 5            | FDA Name/SMILES/SMARTS       | 27           | (29)                     | N/A                       |
| FDA Color Additive List 6            | FDA Name/SMILES/SMARTS       | 29           | (29)                     | N/A                       |
| FDA Color Additive List 7            | FDA Name/SMILES/SMARTS       | 37           | (29)                     | N/A                       |
| Constituents of Cannabis Sativa      | Name/SMILES/SMARTS           | 394          | (30)                     | N/A                       |
| Phytocanniboids                      | Name/SMILES/SMARTS           | 111          | (31)                     | N/A                       |
| Organophosphorous Nerve Toxic Agents | Name/SMILES/SMARTS           | 14           | (32)                     | N/A                       |
| Cengage Bronsted Acids               | Name/SMILES/SMARTS           | 42           | (33)                     | N/A                       |
| Common Regex Patterns                | Mol2                         | 1            |                          | N/A                       |

<p align="center">
  <i>Table 1: GlobalChem Object List Columns: "Chemical List" is the name of the node that contains the chemical list, "Languages" specifies the name and their respective translations, "Number of Entries" is how many molcules exist within one node, "References" are the what resource the molecules were recorded from, and the last column "CGenFF Errors" is the many times CGenFF skipped a molecule. If the value is "N/A" it means it was a node added after testing and allows room for additional chemical space exploraton.</i>
</p>

At the time of writing the list of features include those shown in Table 2. The list range from well defined algorithms implemented into Global-Chem and their respective description and discipline.

| Software Feature              | Description                                                                    |  Code Length | Discipline | Reference |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|---------- |-----------|--------------|
| Validating SMILES             | An adapter to other SMILES platforms (RDKit, PySMILES, SELFIES, PartialSMILES, DeepSMILES, MolVS) to validate by interoperability           | 107        | Cheminformatics          | 
| Decoding Fingerprints and Classifying SMILES            | Decoding your fingerprints to your complex SMILES and to an IUPAC using a annotated dictionary of bit vectors | 129        | Cheminformatics          | 
| SMILES Bidirectional PDF Parsing      | Converting lists of SMILES to 2D drawings in PDF parsing and parsing PDF back to SMILES | 685        | Cheminformatics          | 
| Drug Design Filtering      | Filtering lists of SMILES by a variety of common drug filters (Lipinski Rule of 5, Ghose, Veber, Rule of 3, REOS, Drug-Like, Filters | 137        | Cheminformatics           | 
| Deep Layer Scattering     | Scattering Nodes of Collective SMILES and their relations to each other in Parallel Coordinate Diagram implemented in Plotly | 184        | Cheminformatics          | 
| SMARTS Identifier     | A web application implemented in Flask to test the validation of SMARTS submatching of strings against the MiniFrag Database | 307        | Cheminformatics          | 
| Protonating SMILES    | A distributable version of the Dimorphite-DL package to protonate SMILES over a range of pH with a control over variant production | 56        | Cheminformatics          | 
| Sunbursting SMILES    | Applying a sunburst plot to large collection of SMILES to identify functional groups and pairs of functional groups within the set | 253        | Cheminformatics          | 
| Peptide Sequence to SMILES    | An evolution of Cocktail-Shaker to include Lanthipeptides and covalent sulphur linkages in SMILES strings | 147        | Cheminformatics     
| Visualization SMARTS    | A python application programming interface to port the SMARTS.plus visualizer for SMARTS strings into a jupyter notebook | 47        | Cheminformatics          | 
| One-Hot Encoding SMILES    | A GlobalChem encoder that encodes SMILES for Machine Learning including the '&' denoted as a polymer ex. Diamond | 112        | Cheminformatics          | 
| Principal Component Analysis on SMILES    | A principal component analysis on a list of SMILES with hyperparamter tuning for morgan fingerprinting provided and visualization with Bokeh | 154        | Cheminformatics          | 
| Networkx Adapter    | A graph to graph network adapter between GlobalChem and NetworkX for ease of interoperability for data engineering | 65        | Cheminformatics          | 
| Scaffold Graph Adapter   | An adapter to take a large collection of GlobalChem Nodes and analyze their Structure Hierachy with Scaffold Graphs | 97        | Cheminformatics          | 
| GlobalChem Protein  | An adapter to biopandas to process pdb protein files as well as an implementation of the Bostrom Algorithm to Structurally Filter SMILES | 467        | Bioinformatics          | 
| GlobalChem RNA  | Conversion of RNA Sequence to SMILES and a visualizer for RNA sequences for Python Jupyter Notebooks | 181        | Bioinformatics          | 
| GlobalChem DNA  | Conversion of DNA Sequence to SMILES and a visualizer for DNA sequences for Python Jupyter Notebooks | 181        | Bioinformatics          | 
| GlobalChem Bacteria  | A python model with attributes for general bacteria classifications as well as a common list | 214        | Bioinformatics          | 
| GlobalChem Monoclonal Antibodies  | A python model with attributes for general monoclonoal antibodies classifications as well as a common list | 20        | Bioinformatics          | 
| Z-Matrix Store  | A python model store where users can pull standard z-matrices for molecules queried by their IUPAC | 159        | Quantum Chemistry          | 
| Psi4 Parser  | A python model for analyzing psi4 output files and plotting interaction energy data automatically through Plotly | 193        | Quantum Chemistry          | 
| Moly Adapter  | A software adapter and enhanced functionality for Moly and visualizing HOMO/LUMO orbitals of molecules | 87        | Quantum Chemistry          | 
| GlobalChem Molecule  | A Global Molecule that can parse SMILES, GAFF, CGenFF Stream files into Pandas dataframes, a visualizer with Atom Types and SMILES in RDKit, new mix of cross-discipling languages (SMILES and CGenFF Atom Types) using CXSMILES, CXSMARTS, and CurlySMILES | 386        | Force Fields         | 
| CGenFF Molecule  | A CGenFF Parser that can parse, write edit, and update stream files with Pandas DataFrames | 532        | Force Fields         | 
| GAFF2 Molecule   | A GAFF2 Parser that can parse, write edit, and update stream files with Pandas DataFrames  | 454        | Force Fields         | 
| CGenFF Disimiliarity Score   | A CGenFF dissimilarity algorithm based on the atom types and their tuples of bonded parameters (bonds, angles, dihedrals, impropers) to determine a dissimilarity score | 191        | Force Fields         | 
| Open Source Database Monitor   | An open source database monitor that performs heartbeat checks on common chemical lists running on cloud web servers | 95        | Development Operations         | 
| Plotly Templates   | A Graphing template to use for Ploty to make your data look "pretty" | 80        | Graphing         | 

# Chemical List Selection

Compound lists in Global-Chem can be used to identify specific regions of chemical space that have limited coverage. Therefore, the compound lists in Global-Chem represent future regions of chemical space for force field development. In the CGenFF program we can use larger penalities to indicate a lower extent of analogy to known parameters, information that may be used to identify molecules for additional force field optimization. We passed a variety of Global-Chem objects individually into the software and plotted penalty score distributions of their bonded and non-bonded parameters shown in `Figure 2`. As may be seen the extent of penalties differs significantly for the various lists. Based on the compounds used in the development of CGenFF, we expected the penalties to be lower on molecules that are declared as drugs (Schedule One US Narcotics) and drug-like species (BRAF Kinases Inhibitors for Cancer,  Privileged Scaffolds) whereas we expect the penalty score will be higher for compounds for things that are were not it's original intention ( Emerging PerfluoroAlkyls for Environmental Hazards). 

<p align="center">
  <img width="1000" height="950" src="../images/figures/figure_5.png">
  <i>Figure 2: Penalty Score Probability Distributions</i>
</p>

For force field paramitirization, we want to focus on the most interesting compounds based on human expertise and on the charge
penalty score as an indicator to a new chemical environment that was left unaccounted for in `CGenFF`. This avoids brute force paramitirization 
on mass molecular datasets with no clear intention and allows our force field to remain the best chemical space representation of the community
with human guided direction. To demonstrate our versatility, we use our own tools to explore both avenues of explored and unexplored chemical space.
In Figure 2, when evaluating the distributions of the penalty scores and the nodes that accommodate it, we can begin to evaluate trends to provide
an initial guess of where to look for most likely a compound that wasn't accounted for rapidly. We recorded partial G, V, and A-series toxic agents according 
to Dr. Mirzayanov's account of the Novichok Program (39). Novichok-5 and Sarin contained fluorophosphane bonds that CGenFF has not seen before evident by the partial charge of the phosphorous with penalties upwards of 200. 
We can attribute these high scores to its unique chemical environment specifically designed for warfare which qualifies them as relevant candidates given their history.
If we investigate the covalent warhead inhibitors with a distribution of charge penalty scores ranging from 0 to 300 in Figure 2.
The range indicates that we have accounted for some warheads but in the most recent years they have gained popularity and are being employed elsewhere especially synthesis. Aziridine, in Figure 3, was chosen because it's similar to an epoxide, acting as a good electrophilic drug fragment to extend a molecule by 2 carbons and provide a terminal reactive amine functional group for further extension.
The difference in the CGenFF atom type assignment of the Aziridine is the categorical layer where the oxygen in a 3-ring membered ring system have their own specific sub category and 3-membered ring nitrogen do not. This suggests
the developers of the CGenFF should consider adding a new atom type for this compound. Perfluorobutanoic acid, in Figure 3, is a common herbicide that is popular for the environmental protection agency with a unique chemical pattern where fluorine takes the place of each hydrogen along the alkyl chain, and with, typically, a tail carboxylic acid. This is not typical feature of drug-like molecules and is shown with majority of the charge penalty scores being higher on average on the charge penalty score distribution in Figure 2. This makes them a prime candidate for force field paramitirization. 
Look at the common vitamin list in Figure 2 we found Vitamin C to have a high dihedral penalty scores which was unexpected initially due to it's legacy as an important ubiquitous compound but deeper look into the parameters suggest that the chemical environment,the sp2 carbon in the furanone with surrounding hydroxyls as seen in Figure 3, suggest that this would be tough to navigate through paramitirization. Finally, a dithiolane fragment was discovered using the sunbursting, principal component analysis, and pdf parsing features provided in Global-Chem, when analyzing partial of the Enamine Database, and using Global-Chem's RingsInDrugs as a reference for filtering common ring systems to display emerging ring systems of interest to find **rare** viable ring systems. The 1,3-dithiolane was selected because of its simple cyclopentyl architecture but two sulphurs separated by 1 carbon makes a great potential drug fragment anchor to enhance binding affinity due to the nucleophicility of the sulphurs and their respective geometry. `CGenFF` did not recognize this structure and was ultimately selected as the first compound for force field paramirization with the full account provided in the Supporting Information. 

<p align="center">
<img width="598" alt="Screen Shot 2022-07-19 at 8 09 17 AM" src="https://user-images.githubusercontent.com/11812946/179746963-fc09c731-4f05-4c12-9278-cbccc3ac5efe.png">
  <br>
  <i>Figure 3: MolCloud of Chemical Selection</i>
</p>

# Paramitirization of 1,3 Dithiolane

We truncated dithiolane from the amide and passed through CGenFF (Full data available in the Supporting Information) which indicated that the dilemma was in part due to the extent of puckering caused by the 2 Sulphur atoms within the constrained cyclopentane ring system. T
To begin our parametirization process we chose to focus on `S1-C3-C4-S2`, backbone to the cyclopentane ring and the dihedral from the methyl to one carbon on the backbone `C1-C2-S1-C3`. Since the molecule is symmetric, it makes the complexity of the molecule decrease twofold. 
The parametirization of 1,2-dithiolane was performed using FFParam [Reference Here] following the FFParam Workflow [Reference Here]. To begin our process, we first subject the compound to quantum mechanics (QM) geometry optimization, with Gaussian [Reference HEre], using mp2 theory to treat electron correlation [Reference Here] and basis set of “6-31/+G*” to handle the orbital polarizability of the sulphur atom. 
Our intended goal is to use the QM as a reference target data that the molecular mechanics (MM), CHARMM [Reference Here], should approximately match. We perform potential energy surface (PES) scans around our selected dihedrals and compare the surface of the QM vs the MM. 

To match the PES scan for the MM to the QM we have to tweak “tunable” parameters as defined in charmm potential energy function (i.e force constants, multiplicity) [Reference Here] until we reach a reasonable surface scan and numbers that make common sense. To determine the partial charges, we observe 
the dipole moment induced by the interaction between the atom of interest and water. When the dipole moment of the QM and MM reach within a range 
(< 0.5kcal/mol) we consider that reasonable. 

To accomplish our parametirization we applied the following: for `S1-C3-C4-S2`, if we break the connection ring component 
around the C3-C4 single bond in the  atom ring we obtain a natural rotation of a thiomethyl group. Additional multiplicities of 1 and 2 of varying force constants
seemed to have a negative effect. We added a relatively high force constant of a value of 2.3800 to because this particular 
dihedral is part of a ring where there is a significant energy barrier of rotation due to constraint of the cyclopentane backbone. 

For C1-C2-S1-C3, still maintained the multiplicity of 3 but with a far less reduced force constant of 1.1000. 
This was due to the methyl that replaced the amide allowing some degrees of rotation but the S1 is still constrained within the ring system. 
Final PES scans are displayed in Figure 4. 

<p align="center">
  <img width="600" height="350" src="../images/figures/official_figure_8.png">
  <br>
  <i>Figure 4: Final Potential Energy Scans of dihedrals S1-C3-C4-S2 and C1-C2-S1-C3</i>
</p>

Lastly, the S1-S2 charges needed adjustment. We used Monte Carlo Simulated Annealing (MCSA) method [Reference Here] utilized in
FFparam to predict the approximate partial charges. The sulphur atoms were adjusted to have a partial negative charge of -0.208.

# Theory 

#### GlobalChem Molecule Language

CGenFF and SMILES are built on the same language philosophy yet are independent of each other. Global-Chem serves as a basis generator in combining the languages into something is intuitive to read. CurlySMILES is a subset language of SMILES used to embed a meta data next to a alpha element character for example "C" which means carbon can be read as "C{CG2R61}" a aromatic benzene sp2 carbon. When applying this feature to a more complex molecule we can see how the new bridged language unfolds. We present the first Global-Chem Moleculer Language that contains both CGenFF Atom-Types and SMILES based on scientific inclusion not exclusion (41):

| Molecule                     | Proposed GlobalChem Language                                                                                         | 
|------------------------------|----------------------------------------------------------------------------------------------------------------------|
| Perfluorobutanoic acid       | F{FGA2}C{CG312}(F{FGA2})(C{CG312}(F{FGA2})(C{CG2O2}(O{OG311})=O{OG2D1})F{FGA2})C{CG302}(F{FGA3})(F{FGA3})F{FGA3}     |
| Vitamin C                    | C{CG321}(C{CG311}(C{CG3C51}1C{CG2R51}(=C{CG2R51}(C{CG2R53}(=O{OG2D1})O{OG3C51}1)O{OG311})O{OG311})O{OG311})O{OG311}  |
| Aziridine                    | N{NG311}1C{CG3C31}C{CG3C31}1                                                                                         |
| 1,3-Dithiolane               | C{CG331}C{CG3C51}2S{SG311}C{CG3C52}C{CG3C52}S{SG311}2                                                                |

Using this new language, we can probably determine easily from which atom type could be incorrectly misassigned without looking at the partial charges in conjunction with the SMILES allowing intuition to supersede the penalty score. For example, a N1 in a 3 membered ring is mostly likely not going to be NG311 but probably a new atom type like NG3C31 according to the CGenFF nomenclature. 

#### Decoder Engine

#### Education


# Conclusion

Beyond accessing SMILES strings we've shown the utility of `Global-Chem` to interogate the coverage of the force field `CGenFF`. By partitioning chemical space into well-defined chemical lists, `Global-Chem` allows for regions of chemical space where the CGenFF programs fails or assigns parameters of low analogy to be readily identified. This information will allow for decisions to be made concerning the addition of molecules in the CGenFF training set thereby allowing for systematic improvements in the force field. Accordingly, `Global-Chem` has a number of potential purposes, including teaching and cheminformatics, but our main perogative is to create a free record collection. As `Global-Chem` requires direct user input, if we plant the seed now then, hopefully, our tree will grow. The actual growth of the tree will be decided on by the experts of the  community dedicated to their field. Enjoy.

# Acknowledgements

Thank you to Tyree Wilson, Garrick Centola and Paul Shapiro for their helpful discussions into the usability, functionality, and guidance into the manuscript of Global-Chem. Thank you to Payal Chatterjee for discussions of the 1,3-dithiolane paramirization. Thank you to Steven Fletcher for our discussion on aziridine during a organic synthesis lecture. Appreciation to past mentors James Ryan, Robert Zeigler, and Blake Printy for discussions on good manufacturing practices of python packaging and distribution. Appreciation to the University of Maryland Baltimore, School of Pharmacy, Department of Pharmaceutical Chemistry for promoting a collaborative and useful space for academics from all different scientific disciplines. Financial support from the NIH (GM131710) is acknowledged.

# References

(1) Gorgulla, C.; et. al. An Open-Source Drug Discovery Platform Enables Ultra-Large Virtual Screens. Nature 2020, 580 (7805), 663–668. https://doi.org/10.1038/s41586-020-2117-z.

(2) Irwin, John. J.; et. al. ZINC20—A Free Ultralarge-Scale Chemical Database for Ligand Discovery. ACS Publications 2013, 60 (12), 6065–6073. https://doi.org/10.1021/acs.jcim.0c00675.

(3) Roberts, R. J. PubMed Central: The GenBank of the Published Literature. Proceedings of the National Academy of Sciences of the United States of America 2001, 98 (2), 381–382. https://doi.org/10.1073/pnas.98.2.381.

(4) Cooke-Fox, D. I.; et al. Computer Translation of IUPAC Systematic Organic Chemical Nomenclature. 1. Introduction and Background to a Grammar-Based Approach. Journal of Chemical Information and Computer Sciences 1989, 29 (2), 101–105. https://doi.org/10.1021/ci00062a009.

(5) Weininger, D.; et al. SMILES, a Chemical Language and Information System. 1. Introduction to Methodology and Encoding Rules. Journal of Chemical Information and Computer Sciences 1988, 28 (1), 31–36. https://doi.org/10.1021/ci00062a009.

(6) Fulmer, G. R.; et al. NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist. Organometallics 2010, 29 (9), 2176–2179. https://doi.org/10.1021/om100106e.

(7) Daylight Theory. OpenSmiles.

(8) Lide, D. R.; et al. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data; CRC Press, 2004. https://doi.org/10.5860/choice.37-4225.

(9) Taylor, R. D.; et al. Rings in Drugs. Journal of Medicinal Chemistry 2014, 57 (14), 5845–5859. https://doi.org/10.1021/jm4017625.

(10) Broughton, H. B.; Watson, I. A. Selection of Heterocycles for Drug Design.; PubMed, 2004; Vol. 23, pp 51–58. https://doi.org/10.1016/j.jmgm.2004.03.016.

(11) Welsch, M. E.; et al. Privileged Scaffolds for Library Design and Drug Discovery.; PubMed, 2010; Vol. 14, pp 347–361. https://doi.org/10.1016/j.cbpa.2010.02.018.

(12) Gehringer, M. E.; Laufer, S. A. Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology; ACS Publications, 2019; Vol. 62, pp 5673–5724. https://doi.org/10.1021/acs.jmedchem.8b01153.

(13) Hiorns, R. C. A Brief Guide to Polymer Nomenclature (IUPAC Technical Report).; 2012; Vol. 84, pp 2167–2169. https://doi.org/10.1351/PAC-REP-12-03-05.

(14) Takeuchi, K.; et al. R-Group Replacement Database for Medicinal Chemistry.; 2021; Vol. 7, p FSO742. https://doi.org/10.2144/fsoa-2021-0062.

(15) Petri, L.; et al. An Electrophilic Warhead Library for Mapping the Reactivity and Accessibility of Tractable Cysteines in Protein Kinases.; 2020; Vol. 207, p 112836. https://doi.org/10.1016/j.ejmech.2020.112836.

(16) Hu, H.; et al. Systematic Comparison of Competitive and Allosteric Kinase Inhibitors Reveals Common Structural Characteristics.; 2021; Vol. 214, p 113206. https://doi.org/10.1016/j.ejmech.2021.113206.

(17) Agianian, B.; Gavathiotis, E. Current Insights of BRAF Inhibitors in Cancer.; ACS Publications, 2018; Vol. 61, pp 5775–5793. https://doi.org/10.1021/acs.jmedchem.7b01306.

(18) Isidro-Llobet, A.; et al. Amino Acid-Protecting Groups.; ACS Publications, 2009; Vol. 109, pp 2455–2504. https://doi.org/10.1021/cr800323s.

(19) Pelch, K.; et al. PFAS Health Effects Database: Protocol for a Systematic Evidence Map.; Science Direct, 2019; Vol. 130, p 104851. https://doi.org/10.1016/j.envint.2019.05.045.

(20) Orr, A.; et al. Combining Experimental Isotherms, Minimalistic Simulations, and a Model to Understand and Predict Chemical Adsorption onto Montmorillonite Clays.; PubMed, 2021; Vol. 6, pp 14090–14103. https://doi.org/10.1021/acsomega.1c00481.

(21) ECFR :: 21 CFR Part 1308 - Schedules.

(22) Shulgin, Alexander T., and Ann Shulgin. Pihkal: A Chemical Love Story. 1. ed., 8. print, Transform, 2010.

(23) Vaithianathan, Soundarya, et al. “Effect of Common Excipients on the Oral Drug Absorption of Biopharmaceutics Classification System Class 3 Drugs Cimetidine and Acyclovir.” Journal of Pharmaceutical Sciences, vol. 105, no. 2, Feb. 2016, pp. 996–1005. PubMed

(24) https://github.com/geekan/HowToLiveLonger

(25) Data Abstracted from https://labels.fda.gov/

(26) EXSENS-USA.com. “Lube Lessons 3: The Sex Lube Ingredient Glossary.” EXSENS-USA.Com,

(27) Research, Center for Drug Evaluation and. “Tainted Sexual Enhancement Products.” FDA, June 2022

(28) Belot, Laure  "Alimentation : face aux doutes, les internautes s'organisent". Le Monde.

(29) https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list

(30) Turner, C. E., et al. “Constituents of Cannabis Sativa L. XVII. A Review of the Natural Constituents.” Journal of Natural Products, vol. 43, no. 2, Apr. 1980, pp. 169–234. PubMed

(31) Hanuš, Lumír Ondřej, et al. “Phytocannabinoids: A Unified Critical Inventory.” Natural Product Reports, vol. 33, no. 12, Nov. 2016, pp. 1357–92. PubMed,

(32) Mukherjee, Sudisha, and Rinkoo Devi Gupta. “Organophosphorus Nerve Agents: Types, Toxicity, and Treatments.” Journal of Toxicology, vol. 2020, Sept. 2020, p. 3007984.

(33) PKa Values for Organic and Inorganic Bronsted Acids at 25o Ca.

(34) Mackerell, Alexander D. “Empirical Force Fields for Biological Macromolecules: Overview and Issues.” Journal of Computational Chemistry, vol. 25, no. 13, Oct. 2004, pp. 1584–604. PubMed,

(35) Wang, Junmei, et al. “Development and Testing of a General Amber Force Field.” Journal of Computational Chemistry, vol. 25, no. 9, July 2004, pp. 1157–74. PubMed

(36) Jorgensen, William L., and Julian Tirado-Rives. “The OPLS [Optimized Potentials for Liquid Simulations] Potential Functions for Proteins, Energy Minimizations for Crystals of Cyclic Peptides and Crambin.” Journal of the American Chemical Society, vol. 110, no. 6, Mar. 1988, pp. 1657–66.

(37) Vanommeslaeghe, K., E. et al. CHARMM General Force Field (CGenFF): A Force Field for Drug-like Molecules Compatible with the CHARMM All-Atom Additive Biological Force Fields.” Journal of Computational Chemistry, vol. 31, no. 4, Mar. 2010, pp. 671–90. PubMed Central

(38) Vanommeslaeghe, K., E. Prabhu Raman, et al. “Automation of the CHARMM General Force Field (CGenFF) II: Assignment of Bonded Parameters and Partial Atomic Charges.” Journal of Chemical Information and Modeling, vol. 52, no. 12, Dec. 2012, pp. 3155–68. 

(39) Chemical Weapons Disarmament in Russia: Problems and Prospects. Henry L. Stimson Center, 1995.

(40) Vanommeslaeghe, K., and A. D. MacKerell. “Automation of the CHARMM General Force Field (CGenFF) I: Bond Perception and Atom Typing.” Journal of Chemical Information and Modeling, vol. 52, no. 12, Dec. 2012, pp. 3144–54

(41) Mobley, David L., et al. “Escaping Atom Types in Force Fields Using Direct Chemical Perception.” Journal of Chemical Theory and Computation, vol. 14, no. 11, Nov. 2018, pp. 6076–92.

# Conflict of Interets

ADM is cofounder and CSO and SJ is Commercial Development Director of SilcsBio LLC. Chris Burke is Senior DevOps Engineer at L7 Informatics. Mike Woster is the Chief Revenue Officer of the Linux Foundation, Rebecca Pinette-Dorin is marketing researcher at Exsens.
