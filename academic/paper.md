---
title: 'Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities'

authors:
  - name: Suliman Sharif [first author]
    orcid: 0000-0002-1342-9258
    affiliation: 1
  - name: Ruibin Liu
    orcid: 0000-0001-8395-9353
    affiliation: 1
  - name: Bettina Lier
    affiliation: 4
  - name: Asuka Orr
    orcid: 0000-0003-4628-526X
    affiliation: 1
  - name: Aziza Frank
    affiliation: 1
  - name: Anastasia Croitoru
    orcid: XXX
    affiliation: 5
  - name: Daniel Khavrutskii
    orcid: XXX
    affiliation: 5
  - name: Sunhwan Jo
    affiliation: 2
  - name: Jacob Weiner
    orcid: XXX
    affiliation: 1
  - name: Mingtian Zhao
    orcid: XXX
    affiliation: 1
  - name: Aarion Romany
    affiliation: 1
  - name: Ian Jones
    affiliation: 1
  - name: Shaoqi Zhan
    orcid: 0000-0002-6383-1771
    affiliation: 3
  - name: Chris Burke
    affiliation: 2
  - name: Takayuki Serizawa
    affiliation: 6
  - name: Jared Deacon
    affiliation: 6
  - name: Anmol Kumar
    affiliation: 6
  - name: Mike Woster
    affiliation: 6
  - name: Rebecca Pinette-Dorin
    affiliation: 6
  - name: Elena Yi Chow
    orcid: 0000-0001-5559-1635
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
  - name: University of Vienna
    index: 4
  - name: École Polytechnique
    index: 5
  - name: Takeda Pharmaceuticals
    index: 6
---

# Introduction

The in-silico chemical universe is expanding rapidly as open access titan databases (Enamine Database (20 Billion) [1],
Zinc Database (2 Billion) [2], PubMed Database (68 Million) [3] and cheminformatic tools
to process, manipulate, and derive new compound structures are established. While this chemical data big bang has yielded useful ultra-large datasets they are based on ambiguous classification systems making it difficult to systematically organize them for specific uses.

<p align="center">
  <img width="400" height="300" src="../images/figures/figure_4.png">
  <br />
  <i>Figure 1: Screenshot of the ZincDB request URLS</i>
</p>

For example, in `Figure 1`, the directory setup for downloading ZincDB molecules is shown. As is evident, the information content of the directory nomenclature does not contain information on the compounds they contain, making it nearly impossible to access specific molecules or classes molecules.  Towards overcoming this, partial organizational attempts were made in PubMed, filling chemical data linkages for computational toxicology called Actor for a specific
refactored and refined effort [23]. In another example, for the EnamineDB a scaffold associated with biological activity was designed to target 
Toll-Like Receptors in an object-oriented fashion [24]. However, these organizational methods are difficult
to extend to other systems and can be difficult to implement given the large amount of data.
In addition, the information content of these papers is of limited utility to the common developer. 

To organize chemical compounds we apply the idea of communication. Humans use symbols and drawings to communicate, a set of symbols and the rules to combining them are called a language. Languages can be employed to carry relevant, distinct features and mean something to their respective community, diagrammatically shown in `Figure 2`. 
International Union of Pure and Applied Chemistry (IUPAC) was a coalition that formed in the 1800s and their method of communication was named after the organization, IUPAC. 
IUPAC is a written language that predates even drawing atoms as a method of communication between chemists [6]. 
Other chemical sub-communities adopted the IUPAC language and applied it to their fields that are comprised of different dialects i.e polymer chemistry, organo-metallic chemistry.
Due to it's "first to market" status, the scientific chemical language IUPAC is the legacy language that is the lexical key to unlocking information about a chemical pattern or group. 
But there are problems with the language due to it's length in describing bigger molecules and failure to give a more natural name for smaller molecules ex. water is oxidane. Simply, IUPAC names in organic chemisty papers are impractical, effecting extending the length of a manuscript, while being of limited value given the challenge of interpreting such names. Over time, the IUPAC names naturally turned into a slang ("preferred") language due to humans wanting to speak it while communicating with each other. Effectively, the **natural** chemical language that is extant today is a blend of both a formal and informal nomenclature. 

To compact information, chemists presented drawings of chemical structures but information in such a format is hard to store precisely. Alternatively, SMILES [7] has become a popular 1-D language amongst cheminformaticians as a sufficient way to write and retain 2D chemical connectivity information with ease.  Algorithms
have been designed to abstract and interpolate skeletal patterns and languages from chemical drawings and convert them into SMILES for data processing and analysis. 
A number of these tools, which work to varying degrees of accuracy, have been well summarized by the Blue Obelisk Society Open Source Review [25]. 
Efforts to improve these tools recently have included machine learning (ML) methods that essential "sit" on top of the underlying algoritm to fix any inaccuracies of the method. 
As an alternative we can take another direction, where data is selectively aggregated based on known classifications, popularity and utility, being organized to a degree of functionality that facilitates more widespread use. However, the criteria for such an aggregation of data is built upon human expertise, requiring input from a variety of people to attain the broadness and accessibility that would facilitate scientific discovery. In other words, in the context of a well-classified natural chemical database the major challenge is the enormity of the chemical universe, requring a range of chemical expertise to put together well-thought chemical lists of compounds relevant to their respective communities. Thus, it is necessary to create a tool to allow for a large number of participants to contribute in order for such a data compilation to grow. 
However, most software and especially old software can be difficult to install and handle on top of modern technology thus hindering participation. This situation drives the
need for a tool that is sustainable and readily accessible to potential participants, allowing the database to naturally grow.
This need motivated the development of the presented `Global-Chem` database tool.

To implement `Global-Chem` we selected a coding language that has the ability to write easy objects for particpants to understand; Python [10.5555/159351][7].

<p align="center">
  <img width="700" height="450" src="../images/figures/figure_2.png">
  <br />
  <i>Figure 2: Language organized by category and functionality </i>
</p>

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
  <img width="800" height="700" src="../images/figures/figure_1_official.png">
  <br/>
  <i>Figure 3: Node Network of Global-Chem</i>
</p>

The tree network follows a simple object-oriented pythonic design in conjunction with literature where head nodes are the major corresponding scientific field (example: "Medicinal Chemistry") and their corresponding child nodes are the manuals, articles or books that are the references for the lists.
Each reference object has either the functional groups that correspond to that paper's overall functionality in IUPAC, Preferred Name, Acronyms, SMILES, or SMARTS
format. The motivation for this design was that as more users contribute they can expand into different directories, add their own directory, 
and provide their chemical list of interest. Each paper that is submitted is converted into a `namespace` module, an object
whose name is indicative of it's functionality. An example for the drug design community is the paper "Rings In Drugs" [11] whose
python object equivalent is now "RingsInDrugs" with two functional methods that retrieve the  IUPAC:SMILES/SMARTS dictionary that was embedded included in the master object `Global-Chem`. 
Users can choose to cross reference leaf nodes between each other and do comparative chemical list studies since the IUPAC name and SMILES name are consistent across lists.
Note that not all the SMILES being portrayed are canonical given that users can create their own SMILES, which are not unique. To account for this users can parse `Global-Chem` SMILES into the `RDKit` parser
for canonical SMILES conversion. 

## Data Collection

References and associated compound lists are selected based on the interests of the scientific contributors.  This should include consideration of relevance to the scientific community. To authenticate and validate SMILES strings we employ interoperability tools to other cheminformatic software to verify it's usability.

The IUPAC/SMILES strings may be abstracted in a variety of methods:

- For IUPAC naming we opted for naming things as they were reported in the literature. If no names were available, then we opted to find a natural name to fill the slot.

- For IUPAC naming, we chose to reduce the complexity of the name by opting to remove as much stereochemistry as made sense. 

- For Polymer IUPAC, the site points were omitted from the name and some of the nomenclature adjusted for preferred names
over traditional. For example: 'yl' to mark site points for polymer connections was removed in favor of reduced english complexity. Site points are marked with a virtual atom that can be installed into the SMILES string with the character '*'.

-  For simple molecules one representation of the SMILES can be directly translated using visual 
inspection. This is typically appropriate for compounds at the beginning of a reported list that contain the most common denominator rings. 

- For complex molecules the image can be redrawn in the free version of ChemDraw and then translated into SMILES. 

- For sources where the SMILES are written and the IUPAC is not known the SMILES are translated into ChemDraw and the name retrieved. 
Note that some of the names may be modified based on human inspection in favor of preferred names. 

- In the case of radicals, some SMILES were adjusted to remove the radical chemical feature as they serve as connection points. However in some cases the radical component was maintained, especially in the case of IUPAC blue book common substituents or instellar space where radicals are more unknown and not as well explored.

- SMARTS strings were adapted from the SMILES using RDKit [@Landrum:2019-5]. 

- Stereochemistry which is represented as '@' and '@@' as R and S respectively are removed from the SMILES in order to reduce complexity.

# Data

At the time of writing the list of objects include those shown in Table 1. The list range from well defined classes of chemicals, such as amino acids, to more diverse lists such as Rings in Drugs. In addition, the languages used for each list are given, along with the number entires in the list and the reference.  In addition, the number of times that compounds in each list fail in the CGenFF program, as discussed below, is given.

| Chemical List                       | Languages                    | # of Entries | References               |  CGenFF Errors            |
|-------------------------------------|------------------------------|--------------|--------------------------| --------------------------|
| Amino Acids                         | IUPAC/SMILES/SMARTS          | 20           | Common Knowledge         | 0                         |
| Essential Vitamins                  | Preferred Name/SMILES/SMARTS | 13           | Common Knowledge         | 0                         |
| Common Organic Solvents             | IUPAC/SMILES/SMARTS          | 42           | [8]                      | 3                         |
| Open Smiles                         | IUPAC/SMILES/SMARTS          | 94           | [9]                      | 10                        |
| IUPAC Blue Book (CRC Handbook) 2003 | Preferred Name/SMILES/SMARTS | 333          | [10]                     | 1 (Excluding Radicals)    |
| Rings in Drugs                      | IUPAC/SMILES/SMARTS          | 92           | [11]                     | 0                         |
| Phase 2 Hetereocyclic Rings         | IUPAC/SMILES/SMARTS          | 19           | [12]                     | 0                         |
| Privileged Scaffolds                | IUPAC/SMILES/SMARTS          | 47           | [13]                     | 0                         |
| Common Warheads Covalent Inhibitors | IUPAC/SMILES/SMARTS          | 29           | [14]                     | 4                         |
| Common Polymer Repeating Units      | IUPAC/SMILES/SMARTS          | 78           | [15]                     | 7                         |
| Common R Group Replacements         | IUPAC/SMILES/SMARTS          | 499          | [16]                     | 15                        |
| Electrophillic Warheads for Kinases | Preferred Name/SMILES/SMARTS | 24           | [17]                     | 0                         |
| Privileged Scaffolds for Kinases    | IUPAC/SMILES/SMARTS          | 29           | [18]                     | 0                         |
| BRAF Inhibitors                     | IUPAC/SMILES/SMARTS          | 54           | [19]                     | 5                         |
| Common Amino Acid Protecting Groups | IUPAC/ACRONYM/SMILES/SMARTS  | 346          | [20]                     | 41                        |
| Emerging Perfluoroalkyls            | IUPAC/SMILES/SMARTS          | 27           | [21]                     | 1                         |
| Chemicals For Clay Adsorption       | IUPAC/SMILES/SMARTS          | 33           | [22]                     | 0                         |
| Schedule 1 United States Narcotics  | Preferred Name/SMILES/SMARTS | 240          | [26]                     | 1                         |
| Schedule 2 United States Narcotics  | Preferred Name/SMILES/SMARTS | 60           | [26]                     | 1                         |
| Schedule 3 United States Narcotics  | Preferred Name/SMILES/SMARTS | 22           | [26]                     | 1                         |
| Schedule 4 United States Narcotics  | Preferred Name/SMILES/SMARTS | 77           | [26]                     | 0                         |
| Schedule 5 United States Narcotics  | Preferred Name/SMILES/SMARTS | 8            | [26]                     | 0                         |
| PihKal                              | Preferred Name/SMILES/SMARTS | 179          | [Reference Here]         | 0                         |
| Excipients Cimetidine & Acyclovir   | Preferred Name/SMILES/SMARTS | 14           | [Reference Here]         | 0                         |
| HowToLiveLonger	                    | Preferred Name/SMILES/SMARTS | 4            | [Reference Here]         | 0                         |
| Monoclonal Antibodies               | Preferred Name/SMILES/SMARTS | 19           | [Reference Here]         | 0                         |
| Common Lubricants for Sex Wellness  | Preferred Name/SMILES/SMARTS | 38           | [Reference Here]         | 0                         |
| FDA Tainted Sexual Enhancements     | Preferred Name/SMILES/SMARTS | 4            | [Reference Here]         | 0                         |
| Common Food Salts                   | Preferred Name/SMILES/SMARTS | 14           | [Reference Here]         | 0                         |
| FDA Color Additive List 1           | FDA Name/SMILES/SMARTS       | 12           | [Reference Here]         | 0                         |
| FDA Color Additive List 2           | FDA Name/SMILES/SMARTS       | 15           | [Reference Here]         | 0                         |
| FDA Color Additive List 3           | FDA Name/SMILES/SMARTS       | 16           | [Reference Here]         | 0                         |
| FDA Color Additive List 4           | FDA Name/SMILES/SMARTS       | 39           | [Reference Here]         | 0                         |
| FDA Color Additive List 5           | FDA Name/SMILES/SMARTS       | 27           | [Reference Here]         | 0                         |
| FDA Color Additive List 6           | FDA Name/SMILES/SMARTS       | 29           | [Reference Here]         | 0                         |
| FDA Color Additive List 7           | FDA Name/SMILES/SMARTS       | 37           | [Reference Here]         | 0                         |
| Common Regex Patterns               | Mol2                         | 1            |                          | N/A                       |

<p align="center">
  <i>Table 1: GlobalChem Object List</i>
</p>

# Tests & Applications

A total collection of 2572 unique IUPAC/Preferred Name/Acronym to SMILES/SMARTS was collected (with redundacy) across 37 objects in
an organized fashion by subject. The code was refactored extensively to allow for ease of object addition according to subject and functionality.
`Common Regex Patterns` was omitted from the test because it's not a functional group but rather a substring pattern to extrapolate Tripos `mol2` file information. To test the utility of these lists with other software tests were performed on three open source platforms to determine 
data interoperability. In addition, such tests can indicate cases in which some of the software implemented should be expanded to
include functional groups that could not be parsed. 

### Cheminformatics Test

Global-Chem parsed through seven different tools with majority being successful minus diamond represented with an '&' [Reference Here] and fails with all software including RDKit except the GlobalChem Encoder does account for it. (it is not clear) The percentage of passing is as follows: RDKit 100% [Reference Here], DeepSMILES 99.25% [Reference Here], PartialSMILES 85.7% [Reference Here] , SELFIES 100% [Reference Here], MolVS 98.5% [Reference Here], PySMILES 99.8% [Reference Here]. PartialSMILES proved to be the most robust acceptance/rejection tool in identifying misrepresentations of SMILES. 

| Software        | Number of Failed Compounds | Example of Failed SMILES                                   |
|-----------------|----------------------------|------------------------------------------------------------|
| RDKit           | 0                          |                                                            |
| SELFIES         | 0                          |                                                            |
| Indigo          | 8                          | 'CC([Si](C1=CC=CC=C1)C2=CC=CC=C2)(C)C'                     |
| PySMILES        | 5                          | '[a].[Na+].[K+].[Mg+2].[Ca+2].[Cl-]'                       |
| DeepSMILES      | 8                          | 'OC(C1=CC=CC=C1N2)C2=O', 'C1(N2C(CN=C3)=NN=C2)=C3C=CC=C1'  |    
| MolVS           | 24                         | 'n1ccnc1', 'HF', 'O=N1CCCCC1'                              |
| PartialSMILES   | 337                        | '[CH]C', '[N]=[N+]=[N-]'                                   |

<p align="center">
  <i>Table 2: Intereoperability Results of Common Moleculesfrom Global-Chem against different cheminformatic software</i>
</p>

Indigo's encoder was pretty robust and their software allows for a lot of inteoperabiltiy with different software tools (i.e pdf data parsing of SMILES), when faced with the tert-butyldiphenylsilyl protecting group and the SMILES string with the `Si` is not wrapped in a square brackets that specify an element that doesn't have a complete valence shell. For PySMILES, the inclusion of the '[a]' denoting aromaticity for an "aromatic salt" in the database couldn't be processed. Some other encoders have encoded for an aromaticity keyword as specified in the Daylight Technical Documentation [Reference Here].  DeepSMILES was interesting because it failed on specific functional groups as shown in the example with an oxindole and triazolodiazepines that had complex small branch complexities and moetieies that it didn't foresee existing. MolVS had some interesting results where imidazole (and derivatives) failed probably because it expected for a hydrogen perhaps to be explicity stated due to it's varying protonation states. Hydrofluoric acid was something I was expecting but again the hydrogen actually needed to be enforced with a [H] which is not as intuitive. PartialSMILES proved to be the most robust eluding to SMILES that were partially complete and rejected by their criteria. Failures included a ethyl radical and a azido complex stemming from the interstellarspace molecules. 

### Force Field Test

Access to broad collections of chemical groups will be of interest for development of force fields, *("also known as potential energy functions" not aacurate because a force field is a potential energy function and a collection of parameters. During force field development we mostrly work on adding new parameters or ameliorating the existing ones. If the potential energy function is modified, a completely new molecular mechanics model will be created and it would be incompatible with the previous parameters)*,[MacKerell:2004-10] for molecular modeling and molecular dynamic simulations, 
allowing for studies on a wider range of chemicals and biological systems. The ability of a force field to treat molecules in the database can also serve as dual interoperable test 
for SMILES strings. Popular force fields such as General Amber ForceField (GAFF) [Wang:2004-7], Optimized Potentials for Liquid Simulations (OPLS)
[Jorgensen:1988-7], and Charmm General Force Field (CGenFF) [Vanommeslaeghe:2010-3] are based on collections of chemicals that are representative of the particular region 
of chemical space that the force field was designed to cover. In practice, this involves the atom-typing engine of each force field being applied 
to each molecule followed by assignment of the appropriate parameters by analogy to molecules present in the force field. The accuracy of new parameters is limited by the coverage of a specific force field.
Thus, the compound lists in Global-Chem can be used to identify specific regions of chemical space that have limited coverage. Therefore, the compound lists in Global-Chem represent future regions of chemical space for force field development. In the present study, we used CGenFF to check its tolerance level for
the range of molecules currently in Global-Chem. To facilitate this an in-house extension of CGenFF was used that can assign atom types from `SDF` bond type column.
This enabled us to pass the SMILES strings through `RDKit` and transform `SDF` to a `CGenFF` stream output. The resulting failures are also presented in Table 1.
It should be noted by nature of the data processing workflow anything that fails in `RDKit` fails in `CGenFF`.


The CGenFF acted as the foundation for the development of the `CGenFF Program` [Vanommeslaeghe:2010-3] that inputs
molecules and outputs the topology information and parameters required to perform various types of molecular modeling and simulations
using programs such as CHARMM, NAMD, OpenMM and Gromacs [Jo:2008-06]. To test the ability of the CGenFF program to handle the chemical 
lists in Global-Chem each list was individually submitted to the program. As shown in Table 1, the range of failures varies widely.
The majority of lists associated with biological or drug-like molecules have zero failures. In contrast, lists such as Common R-group replacements 
or Protecting Groups show a number of failures. In addition, CGenFF does not cover radicals, which were excluded from the analysis. 
Thus, Global-Chem allows for areas of poor coverage of CGenFF to be identified, information that can be used to facilitate future force field development.

More granular information on the regions of chemical space that need additional development in CGenFF can be made based on the CGenFF penalty score distribution [Vanommeslaeghe:2012].
Penalty scores are attributed to molecules by the CGenFF program whose entire chemical connectivity is not present in CGenFF. When an arbitrary molecule is 
passed through the CGenFF program it navigates through a set of rules that represent a atom type similarity network tree.
Once atom types along with chemical connectivity are known, bonded parameters available in CGenFF are assigned to the molecule. 
If exact matches of the bonded parameters are not available, a second tree traversal browses for alternate parameter by using a second rules files that assigns penalties based on the analogy to known parameters. 
Once the lowest penalty score bonded parameter substitutions are determined, the `CGenFF Program` assigns those parameters along with the associated penalties.
In addition, the program identifies the original parameters that is also output into the stream file  used in the 
various molecular modeling programs. Partial atomic charges and associated penalties are assigned through an extended bond-charge increment scheme. It consist in associating atom type along with chemical connectivity (including  bond, angle and dihedral) with charge increment values subtracted from the atoms formal charge. Thus, while the CGenFF program can successfully ingest a large number of molecules, the majority of those molecules are assigned penalties that indicate the level of analogy of the assigned bonded parameters and charges. Larger penalities indicate a lower extent of analogy to known parameters, information that may be used to identify molecules for additional force field optimization.

Motivated by the availability of the CGenFF penalty scores we passed a variety of objects individually into the `CGenFF program` using our in-house version that can process SMILES strings and recorded the results.
The penalty score distributions are shown in `Figure 6` in a rug fashion using Plotly [Plotly] to show the extent of `CGenFF` penalites
for the different chemical lists. As may be seen the extent of penalties differs significantly for the various lists. 
To understand the utility of this information we focus on five leaf nodes: Schedule One US Narcotics (240), BRAF Kinases Inhibitors for Cancer (54), Privileged Scaffolds (47), Common Warheads (29), [Gehringer:2019-6] and Emerging PerfluoroAlkyls (27). Schedule One are active drugs that are popular in the black market [21CFRPart1], kinase inhibitors should contain drug-like features, privileged scaffolds are selected compounds produced by nature, warheads are designed for covalent inhibition, and PerfluoroAlkyls include herbicides and other compounds that are toxic to humans. Based on the compounds used in the development of CGenFF,
we expected the penalties to be lower on drugs and drug-like species and higher for compounds from chemical manufacturing. 

<p align="center">
  <img width="1000" height="950" src="../images/figures/figure_5.png">
  <i>Figure 4: Penalty Score Probability Distributions</i>
</p>

From `Figure 4`, if we use the charge penalty score as a metric for performance, it is evident that the `CGenFF program` 
assigns parameters with generally low penalty scores, less than 200, for Schedule One and BRAF Kinase Inhibitors owed to its
initial training set of "drug-like" molecules. Privileged Scaffolds encompass a lot of natural products which 
have functional groups that fall into the definition of "drug-like" but not all as indicated by the purple lines  
between penalties 200 and 400 representing high charge penalties. A similar trend is seen with the Common Warheads, with most charge penalties being less than 200, but two prominent purple lines between 200 and 400 associated with high charge penalties, as these compounds contain drug-like features along with reactive functional groups that were not in the CGenFF training set. With both of these lists, it would be useful to identify specific molecules with high penalties and include them in the CGenFF training set. And lastly, Perfluoroalkyls are used in chemical manufacturing of everyday goods [Pelch:2019-9]. While the `CGenFF` training set did include halogens [Soteras:2016-10], motivated by their inclusion in many drugs, `CGenFF` was not extended to perfluoroalkyls.

Accordingly, for this list, there are no low penalty scores with the scores clustered in the intermediate range. This is consistent with halogens being inlcuded the training of CGenFF but the specific connectivity of perfluoroalkyls (long haloalkyl chains) not being included.
Accordingly, if even a few perfluoroalkyls are added to the `CGenFF` training set it will help reduce penalties and improve that treatment of this class of molecules making CGenFF of more utility to the chemical hazard community. 

In addition to the ability of CGenFF to treat the selected chemical lists discussed above other noteworthy failures are explained. For example, cyclobutadiene is a non-traditional ring system with a lot of ring strain although the carbon atom types are common.
`CGenFF` might determine that this particular ring system with it's existing atom type network is not allowed or detrimental to the network if added and needs to be handled with care. An interesting group that fails in CGenFF are allene-based compounds and perhaps warrants extension of the force.
Silicon has not been included in CGenFF leading to the failures of the silicon-based compounds. Similarly, the IUPAC blue book valuable list includes radicals, which are relevant for synthesis purposes. This is another class for `CGenFF`has not yet been parametirized. Full logs of failed compounds are found in the `tests` directory in the github repository. 

# Chemical Selection 

For forcefield parametrization, there are two avenues for sufficient chemical selection. First, **common** compounds relevant to the community which expands the forcefield coverage into relevant chemical space while avoiding compounds that would most likely never exist, this avoids performing brute force parametirization on mass molecular datasets. Second, is **rare** compounds that were valuable history but has been potentially buried in data that we have forgotten about them. To demonstrate the versaility of our software, we will use it to explore both avenues of explored and unexplored chemical space. Let's revisit Figure 4, when evaluating the distributions of the penalty scores and the nodes that accodomate it, we can intuitively guess where to look for most likely a compound that we didn't account for. If we look into the covalent warhead inhibitors with charge penalty scores ranging from 0 to 300. Warhead inhibtors are usually small esoteric chemical environments that are most likely an atom type that CGenFF program hasn't seen before and is most likely misreprensting it. Aziridine, as labeled in Figure 6, is most likely a good candidate because it's similar to epoxide but the atom type assignment is different where oxygens in a 3-ring membered ring system have their own specific sub category and nitrogens do not. Aziridine is also useful in synthtesis as a a great electrophile drug fragment to add 2 carbones and a terminal amine to a molecule. Aziridine has a recent popularity as well with the rise of covalent warhead inbitors  being subject to act by the cysteines on proteins. This makes aziridine a prime candidate for forcefield paramtirization relevant to the chemical community. By applying this analysis method of the charge distribution and casual inference to determine relevance we selected our list.

<p align="center">
  <img width="500" height="350" src="../images/figures/official_figure_6.png">
  <br>
  <i>Figure 5: Chemical Selection List determined off relevance and charge distribution score and ease of paramitirization with their respective atom-types</i>
</p>

To detect rare based compounds we utilize the sunbursting software method within the Global-Chem Extensions package to illunminate a portion of the Enamine Database that passed CGenFF to find new ring system compounds. 
The first layer is a lexical key to the functional groups abstracted from the the "RingsInDrugs" node reported in GlobalChem that could be a good reference to common rings with possible rare atom type scenarios. 
The next layer identifies chemical compounds that have at least two matches with the ring functional groups being considered meaning they exist within the same SMILES string. 
The next layer is the total atom charge penalty score for the SMILES string and the final layer consists of the individual, highest penalty atom and its atom type existing within that SMILES string. 
It is the highest penalized atom type in the functional group that is the information needed by the user to initiate the parameter optimization process. By using this combination of common ring systems with rare atom types we can possibly find new rings that are similar but slightly different according to the atom type.
Which could be useful to the community for small functional group conversion and preserveing functionality while possibly expanding it's application.

<p align="center">
  <img width="900" height="1350" src="../images/figures/official_figure_7.png">
  <br>
  <i>Figure 6: Application of the Sunbursting Method Applied on CGenFF Atom Types and RingsInDrugs GlobalChem Node</i>
</p>

We found an amide in conjunction with the disulphur embedded cyclopentyl ring called officially a dithiolane. 
Dithiolanes are easy to synthesize [Reference Here], naturally occurring [Reference Here], and becoming an emerging potent anchor fragments in drug design [Reference Here] and
 their parameters were deemed the most valuable to append to CGenFF out of the dataset.

# CGenFF ForceField Parametirization of the 1,2-Dithiolane

We truncated dithiolane from the amide and passed through CGenFF (Full data available in the Supporting Information) which indicated that the dilemma was in part due to the extent of puckering caused by the 2 Sulphur atoms within the constrained cyclopentane ring system. T
To begin our parametirization process we chose to focus on `S1-C3-C4-S2`, backbone to the cyclopentane ring and the dihedral from the methyl to one carbon on the backbone `C1-C2-S1-C3`. Since the molecule is symmetric, it makes the complexity of the molecule decrease twofold. 
The parametirization of 1,2-dithiolane was performed using FFParam [Reference Here] following the FFParam Workflow [Reference Here]. To begin our process, we first subject the compound to quantum mechanics (QM) geometry optimization, with Gaussian [Reference HEre], using mp2 theory to treat electron correlation [Reference Here] and basis set of “6-31/+G*” to handle the orbital polarizability of the sulphur atom. 
Our intended goal is to use the QM as a reference target data that the molecular mechanics (MM), CHARMM [Reference Here], should approximately match. We perform potential energy surface (PES) scans around our selected dihedrals and compare the surface of the QM vs the MM. 

To match the PES scan for the MM to the QM we have to tweak “tunable” parameters as defined in charmm potential energy function (i.e force constants, multiplicity) [Reference Here] 
until we reach a reasonable surface scan and numbers that make common sense. To determine the partial charges, we observe 
the dipole moment induced by the interaction between the atom of interest and water. When the dipole moment of the QM and MM reach within a range 
(< 0.5kcal/mol) we consider that reasonable. 

To accomplish our parametirization we applied the following: for `S1-C3-C4-S2`, if we break the connection ring component 
around the C3-C4 single bond in the  atom ring we obtain a natural rotation of a thiomethyl group. Additional multiplicities of 1 and 2 of varying force constants
seemed to have a negative effect. We added a relatively high force constant of a value of 2.3800 to because this particular 
dihedral is part of a ring where there is a significant energy barrier of rotation due to constraint of the cyclopentane backbone. 

For C1-C2-S1-C3, still maintained the multiplicity of 3 but with a far less reduced force constant of 1.1000. 
This was due to the methyl that replaced the amide allowing some degrees of rotation but the S1 is still constrained within the ring system. 
Final PES scans are displayed in Figure 8. 

<p align="center">
  <img width="700" height="450" src="../images/figures/official_figure_8.png">
  <br>
  <i>Figure 7: Final Potential Energy Scans of dihedrals S1-C3-C4-S2 and C1-C2-S1-C3</i>
</p>

Lastly, the S1-S2 charges needed adjustment. We used Monte Carlo Simulated Annealing (MCSA) method [Reference Here] utilized in
FFparam to predict the approximate partial charges. The sulphur atoms were adjusted to have a partial negative charge of -0.208.
The initial and final modeled result of the dithiolane is displayed in Figure 7.

# Conclusion 

`Global-Chem` was developed to facilitate accessing lists of known chemical compounds as objects to allow them to be used in the context of python-based workflows.
However, it can also facilitate the evaluation of other tools to access chemical information in the form of SMILES. An interesting observation from the present data is the ability of tools to handle the ampersand `&` operator in SMILES for materials. 
For example, diamond is a common carbon substance whose SMILES strings is indicated in the OpenSMILES
documentation as a `C&1&1&1&1`. As shown in Table 2, this fails in both `RDKit` and `Indigo` indicating that improved handling of the `&` operator is required. 

Beyond accessing SMILES stings we've shown the utility of `Global-Chem` to interogate the coverage of the force field `CGenFF`. By partitioning chemical space into well-defined chemical lists, `Global-Chem` allows for regions of chemical space where the CGenFF programs fails or assigns parameters of low analogy to be readily identified. This information will allow for decisions to be made concerning the addition of molecules in the CGenFF training set thereby allowing for systematic improvements in the force field.

# Statement of Purpose

`Global-Chem` was developed to facilitate the ability of scientists in both academia and industry to make their compounds of interest readily available to the scientific community in the form of objects that may be directly accessed from python. 
Accordingly, `Global-Chem` has a number of potential purposes, including teaching and cheminformatics, but our main perogative is to create a free record collection.
As `Global-Chem` requires direct user input, if we plant the seed now then, hopefully, our tree will grow. 
The actual growth of the tree will be decided on by the common chemical community and experts in the field. Enjoy. 

# Acknowledgements

Thank you to Tyree Wilson, and Paul Shapiro for their helpful discussions into the usability and functionality of Global-Chem.
Appreciation to past mentors James Ryan, Robert Zeigler, and Blake Printy for discussions on good manufacturing practices of python packaging and distribution.
Appreciation to the University of Maryland Baltimore, School of Pharmacy, Department of Pharmaceutical Chemistry for promoting a collaborative and useful space for academics. Financial support from the NIH (GM131710) is acknowledged.

# References

(1) Gorgulla, C.; et. al. An Open-Source Drug Discovery Platform Enables Ultra-Large Virtual Screens. Nature 2020, 580 (7805), 663–668. https://doi.org/10.1038/s41586-020-2117-z.

(2) Irwin, John. J.; et. al. ZINC20—A Free Ultralarge-Scale Chemical Database for Ligand Discovery. ACS Publications 2013, 60 (12), 6065–6073. https://doi.org/10.1021/acs.jcim.0c00675.

(3) Roberts, R. J. PubMed Central: The GenBank of the Published Literature. Proceedings of the National Academy of Sciences of the United States of America 2001, 98 (2), 381–382. https://doi.org/10.1073/pnas.98.2.381.

(4) Landrum, G. RDKit: Open-Source Cheminformatics Software. 2006. https://doi.org/10.5281/zenodo.591637.

(5) Landrum, G. Indigo. 2009.

(6) Cooke-Fox, D. I.; et al. Computer Translation of IUPAC Systematic Organic Chemical Nomenclature. 1. Introduction and Background to a Grammar-Based Approach. Journal of Chemical Information and Computer Sciences 1989, 29 (2), 101–105. https://doi.org/10.1021/ci00062a009.

(7) Weininger, D.; et al. SMILES, a Chemical Language and Information System. 1. Introduction to Methodology and Encoding Rules. Journal of Chemical Information and Computer Sciences 1988, 28 (1), 31–36. https://doi.org/10.1021/ci00062a009.

(8) Fulmer, G. R.; et al. NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist. Organometallics 2010, 29 (9), 2176–2179. https://doi.org/10.1021/om100106e.

(9) Daylight Theory. OpenSmiles.

(10) Lide, D. R.; et al. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data; CRC Press, 2004. https://doi.org/10.5860/choice.37-4225.

(11) Taylor, R. D.; et al. Rings in Drugs. Journal of Medicinal Chemistry 2014, 57 (14), 5845–5859. https://doi.org/10.1021/jm4017625.

(12) Broughton, H. B.; Watson, I. A. Selection of Heterocycles for Drug Design.; PubMed, 2004; Vol. 23, pp 51–58. https://doi.org/10.1016/j.jmgm.2004.03.016.

(13) Welsch, M. E.; et al. Privileged Scaffolds for Library Design and Drug Discovery.; PubMed, 2010; Vol. 14, pp 347–361. https://doi.org/10.1016/j.cbpa.2010.02.018.

(14) Gehringer, M. E.; Laufer, S. A. Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology; ACS Publications, 2019; Vol. 62, pp 5673–5724. https://doi.org/10.1021/acs.jmedchem.8b01153.

(15) Hiorns, R. C. A Brief Guide to Polymer Nomenclature (IUPAC Technical Report).; 2012; Vol. 84, pp 2167–2169. https://doi.org/10.1351/PAC-REP-12-03-05.

(16) Takeuchi, K.; et al. R-Group Replacement Database for Medicinal Chemistry.; 2021; Vol. 7, p FSO742. https://doi.org/10.2144/fsoa-2021-0062.

(17) Petri, L.; et al. An Electrophilic Warhead Library for Mapping the Reactivity and Accessibility of Tractable Cysteines in Protein Kinases.; 2020; Vol. 207, p 112836. https://doi.org/10.1016/j.ejmech.2020.112836.

(18) Hu, H.; et al. Systematic Comparison of Competitive and Allosteric Kinase Inhibitors Reveals Common Structural Characteristics.; 2021; Vol. 214, p 113206. https://doi.org/10.1016/j.ejmech.2021.113206.

(19) Agianian, B.; Gavathiotis, E. Current Insights of BRAF Inhibitors in Cancer.; ACS Publications, 2018; Vol. 61, pp 5775–5793. https://doi.org/10.1021/acs.jmedchem.7b01306.

(20) Isidro-Llobet, A.; et al. Amino Acid-Protecting Groups.; ACS Publications, 2009; Vol. 109, pp 2455–2504. https://doi.org/10.1021/cr800323s.

(21) Pelch, K.; et al. PFAS Health Effects Database: Protocol for a Systematic Evidence Map.; Science Direct, 2019; Vol. 130, p 104851. https://doi.org/10.1016/j.envint.2019.05.045.

(22) Orr, A.; et al. Combining Experimental Isotherms, Minimalistic Simulations, and a Model to Understand and Predict Chemical Adsorption onto Montmorillonite Clays.; PubMed, 2021; Vol. 6, pp 14090–14103. https://doi.org/10.1021/acsomega.1c00481.

(23) Judson, R. S.; et al. Aggregating Data for Computational Toxicology Applications: The U.S. Environmental Protection Agency (EPA) Aggregated Computational Toxicology Resource (ACToR) System.; PubMed Central, 2012; Vol. 13, pp 1805–1831. https://doi.org/10.3390/ijms13021805.

(24) Perez-Regidor, L.; et al. Virtual Screening Approaches towards the Discovery of Toll-Like Receptor Modulators.; PubMed Central, 2016; Vol. 17, p 1508. https://doi.org/10.3390/ijms17091508.

(25) O-Boyle, N.; et al. Open Data, Open Source and Open Standards in Chemistry: The Blue Obelisk Five Years On.; BioMed Central, 2011; Vol. 3, p 37. https://doi.org/10.1186/1758-2946-3-37.

(26) ECFR :: 21 CFR Part 1308 - Schedules.

(27) Van Rossum, G.; Drake, F. L. Python 3 Reference Manual; CreateSpace: Scotts Valley, CA, 2009.

# Conflict of Interets

ADM is cofounder and CSO and SJ is Commercial Development Director of SilcsBio LLC. Chris Burke is Senior DevOps Engineer at L7 Informatics. 
