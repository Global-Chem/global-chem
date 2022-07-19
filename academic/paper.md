---
title: 'Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities'

authors:
  - name: Suliman Sharif [first author]
    orcid: 0000-0002-1342-9258
    affiliation: 1
  - name: Ruibin Liu
    orcid: 0000-0001-8395-9353
    affiliation: 1
  - name: Asuka Orr
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
    orcid: XXX
    affiliation: 1
  - name: Nathaniel McClean
  - orcid: XXX
  - affiliation: 1
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
to process, manipulate, and derive new compound structures are established. While this chemical data big bang has yielded useful ultra-large datasets they are based on ambiguous classification systems making it difficult to systematically organize them and select the most relevant compounds for forcefield paramirization. These organizational methods are hard to extend to other systems and can be complicated to implement given the large amount of data. In addition, the information content of these papers is of limited utility to the common developer. 

To organize chemical compounds we apply the idea of communication. International Union of Pure and Applied Chemistry (IUPAC) is a written language that predates even drawing atoms as a method of communication between chemists (4). Over time, the IUPAC names naturally turned into a slang ("preferred") language due to humans wanting to speak it while communicating with each other. Effectively, the **natural** chemical language that is extant today is a blend of both a formal and informal nomenclature. To compact information, chemists presented drawings of chemical structures but information in such a format is hard to store precisely. Alternatively, SMILES (5) has become a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D chemical connectivity information with ease. We developed an in-house automated workflow to process SMILES to Charmm General Force Field (CGenFF) Atom Types (40) to enable us to safely exchange chemical information between collaborators by something we can naturally speak to select the most relevant compounds. Algorithms have been designed to abstract and interpolate skeletal patterns and languages from chemical drawings and convert them into SMILES for data processing and analysis, however, we chose to write ours **manually** to effectively learn the language, organize the chemicals appropiately by their purpose, and develop rules to ensure data integrity and consistency. 

Selecting chemical compounds requires expertise. Expertise is gained by experience and reading about a dedicated discipline. Dedicated displines most often has a set of common functional groups that are relevant to that community, this allows us to focus on compounds that are valuable. We do not need all the compounds, since a lot of them are not useful or not possible. In our paper, we describe how we elected compounds and organized them into a open source knowledge graph toolkit, Global-Chem, to allow ease of communication across the community to select the best possible molecule candidates for force field paramitirization. 

<p align="center">
<img width="1000" alt="Screen Shot 2022-07-16 at 5 29 41 PM" src="https://user-images.githubusercontent.com/11812946/179372511-61758864-6b0a-410e-b15f-578fd8227a14.png">
    <br>
    <i>Figure 1: Knowledge Graph of GlobalChem</i>
</p>

# Data

At the time of writing the list of objects include those shown in Table 1. The list range from well defined classes of chemicals, such as amino acids, to more diverse lists such as Rings in Drugs. In addition, the languages used for each list are given, along with the number entires in the list and the reference. The number of times that compounds in each list fail in our automated SMILES to CGenFF workflow is given where if the value is "N/A" it means it was a node added after testing and allows room for additional chemical space exploraton.


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
  <i>Table 1: GlobalChem Object List</i>
</p>

# Chemical Selection for Force Fields

Access to broad collections of chemical groups will be of interest for development of force fields (34), for molecular modeling and molecular dynamic simulations, allowing for studies on a wider range of chemicals and biological systems. The ability of a force field to treat molecules in the graph can also serve as a performance test coverage of useful chemical space general popular force fields such as General Amber Force Field (35), Optimized Potentials for Liquid Simulations (36), and Charmm General Force Field (37) for which they are designed to cover. In practice, this involves the atom-typing engine of each force field being applied to each molecule followed by assignment of the appropriate parameters by analogy to molecules present in the force field. The accuracy of new parameters is limited by the coverage of a specific force field. Thus, the compound lists in Global-Chem can be used to identify specific regions of chemical space that have limited coverage. Therefore, the compound lists in Global-Chem represent future regions of chemical space for force field development. In the present study, we used CGenFF to check its tolerance level for a range of molecules currently in Global-Chem. More granular information of the chemical space regions that need additional development in CGenFF can be made based on the CGenFF penalty score distribution (38). Penalty scores are attributed  by the CGenFF program to molecules whose entire chemical connectivity is not present in CGenFF. When a molecule is passed to the CGenFF program it navigates through a set of rules that represent a atom type similarity network, a tree structure, similar to the Global-Chem file structure. Once its atom types along with chemical connectivity are known, bonded parameters available in CGenFF are assigned to the molecule. If exact matches of the bonded parameters are not available, a second tree traversal browses for alternate parameter by using a second-rules file that assigns penalties based on the analogy to known parameters. Once the bonded parameter substitutions with the lowest penalty score  are determined, the CGenFF Program assigns those parameters along with the associated penalties. In addition, the program identifies the original parameters that are also output into the stream file used in the various molecular modeling programs. Partial atomic charges and associated penalties are assigned through an extended bond-charge increment scheme. It consists in associating atom type along with chemical connectivity (including  bond, angle and dihedral) with charge increment values subtracted from the atoms formal charge. Thus, while the CGenFF program can successfully ingest a large number of molecules, the majority of those molecules are assigned penalties that indicate the level of analogy of the assigned bonded parameters and charges. Larger penalities indicate a lower extent of analogy to known parameters, information that may be used to identify molecules for additional force field optimization.

Motivated by the availability of the CGenFF penalty scores we passed a variety of objects individually into the software and plotted penalty score distributions of their bonded and non-bonded parameters shown in `Figure 2`. As may be seen the extent of penalties differs significantly for the various lists. Based on the compounds used in the development of CGenFF, we expected the penalties to be lower on molecules that are declared as drugs (Schedule One US Narcotics) and drug-like species (BRAF Kinases Inhibitors for Cancer,  Privileged Scaffolds) whereas we expect the penalty score will be higher for compounds for things that are were not it's original intention ( Emerging PerfluoroAlkyls for Environmental Hazards). 

<p align="center">
  <img width="1000" height="950" src="../images/figures/figure_5.png">
  <i>Figure 2: Penalty Score Probability Distributions</i>
</p>

For force field paramitirization, we want to focus on the most interesting compounds based on human expertise and on the charge
penalty score as an indicator to a new chemical environment that was left unaccounted for in `CGenFF`. This avoids brute force paramitirization 
on mass molecular datasets with no clear intention and allows our force field to remain the best chemical space representation of the community
with human guided direction. To demonstrate our versatility, we use our own tools to explore both avenues of explored and unexplored chemical space.
Let's revisit Figure 2, when evaluating the distributions of the penalty scores and the nodes that accommodate it, we can begin to evaluate trends to provide
an initial guess of where to look for most likely a compound that wasn't accounted for rapidly.  We recorded partial G, V, and A-series toxic agents according 
to Dr. Mirzayanov's account of the Novichok Program (39). Novichok-5 and Sarin contained fluorophosphane bonds that CGenFF has not seen before evident by the partial charge of the phosphorous with penalties upwards of 200, 
which we can owe to its unique chemical environment specifically designed for warfare which qualifies them as relevant candidates given their history.
If we investigate the covalent warhead inhibitors with a distribution of charge penalty scores ranging from 0 to 300.
The range indicates that we have accounted for some warheards but in the most recent years they have gained popularity and are being employed elsewhere especiall synthesis. Aziridine was chosen because it's similar to an epoxide 
because of it's usefulness in synthesis acting as a great electrophile drug fragment to add expand an alkyl chain by 2 carbons and provide a terminal amine functional group for chemical extension.
The difference in the CGenFF atom type assignment of the Aziridine is the categorical layer where the oxygen in a 3-ring membered ring system have their own specific sub category and 3-membered ring nitrogen do not. This suggests
the developers of the CGenFF should consider adding a new atom type for this compound. Common herbicides are very popular for the environmental protection agency
with a unique chemical pattern of long alkyl chains surrounded in fluorine halogens with, typically, a tail carboxylic acid. The smallest of the Perfluoroalkanes was perfluorobutanoic acid, where each carbon substituent, when available, is a fluorine. 
This is not typical feature of drug-like molecules and is shown with majority of the charge penalty scores being averagely high on the distribution. This makes them a prime candidate for force field paramitirization. 
What was interesting that the common vitamin list we found Vitamin C to have a high charge penalty score which was unexpected initially but looking more at the structure with the sp2 carbon enwrapped in the furanone might be difficult to add but makes a formidable opponent to add
due to it's legacy as an important ubiquitous compound in human history. Finally, a dithiolane fragment was discovered, when analyzing partial of the Enamine Database, and using Global-Chem's RingsInDrugs as a reference for filtering common ring systems and display new ring systems of interest. 
The 1,3-dithiolane was selected because of it's simple cyclopentyl architecture but two sulphurs separated by 1 carbon but makes a great potential drug 
fragment anchor to enhance binding affinity due to the nucleophicility of the sulphurs and their respective geometry. `CGenFF` did not recognize this structure 
and was ultimately selected as the first compound for force field paramirization with the full account provided in the Supporting Information. 

<p align="center">
<img width="598" alt="Screen Shot 2022-07-19 at 8 09 17 AM" src="https://user-images.githubusercontent.com/11812946/179746963-fc09c731-4f05-4c12-9278-cbccc3ac5efe.png">
  <br>
  <i>Figure 3: MolCloud of Chemical Selection</i>
</p>

# Conclusion 

`Global-Chem` was developed to facilitate the ability of scientists in both academia and industry to make their compounds of interest readily available to the scientific community in the form of objects that may be directly accessed from python. However, it can also facilitate the evaluation of other tools to access chemical information in the form of SMILES. Beyond accessing SMILES stings we've shown the utility of `Global-Chem` to interogate the coverage of the force field `CGenFF`. By partitioning chemical space into well-defined chemical lists, `Global-Chem` allows for regions of chemical space where the CGenFF programs fails or assigns parameters of low analogy to be readily identified. This information will allow for decisions to be made concerning the addition of molecules in the CGenFF training set thereby allowing for systematic improvements in the force field. Accordingly, `Global-Chem` has a number of potential purposes, including teaching and cheminformatics, but our main perogative is to create a free record collection. As `Global-Chem` requires direct user input, if we plant the seed now then, hopefully, our tree will grow. The actual growth of the tree will be decided on by the experts of the  community dedicated to their field. Enjoy.

# Acknowledgements

Thank you to Tyree Wilson, and Paul Shapiro for their helpful discussions into the usability and functionality of Global-Chem. Thank you to Steven Fletcher for our discussion on aziridine during a organic synthesis lecture. Appreciation to past mentors James Ryan, Robert Zeigler, and Blake Printy for discussions on good manufacturing practices of python packaging and distribution. Appreciation to the University of Maryland Baltimore, School of Pharmacy, Department of Pharmaceutical Chemistry for promoting a collaborative and useful space for academics from all different scientific disciplines. Financial support from the NIH (GM131710) is acknowledged.

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

# Conflict of Interets

ADM is cofounder and CSO and SJ is Commercial Development Director of SilcsBio LLC. Chris Burke is Senior DevOps Engineer at L7 Informatics. Mike Woster is the Chief Revenue Officer of the Linux Foundation, Rebecca Pinette-Dorin is marketing researcher at Exsens.
