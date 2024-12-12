<h1 align="center">Organizing The Chemical Universe</h1>

Global Chem is a public dictionary of common chemical lists using the Common Chemical Name as input and SMILES/SMARTS as output organized by their respective community in a knowledge graph. Our hope is this repository serves as a base for the population to govern how the chemicals we use in things like Food, Clothing, Environment, Materials, Drugs, War and a lot more are beneficial for all of us.

For more documentation, please navigate to our Github Wiki.

- Demo https://colab.research.google.com/drive/1HaDvAYpaX0_2Eqx1NSl5uJTRkE_0na-L?usp=sharing
- Join Our Community: https://discord.gg/xygRVwNZP8

![Global Chem Downloads](https://pepy.tech/badge/global-chem)
![Global Chem Extensions Downloads](https://pepy.tech/badge/global-chem-extensions)
![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)

<p align="center">
<img width="800" alt="Screen Shot 2022-07-16 at 5 29 41 PM" src="https://user-images.githubusercontent.com/11812946/179372564-c286b115-af14-4ad8-a37f-0a216297b6c1.png">
</p>

QuickStart
==========

GlobalChem going to be distribute via PyPi as saperate modules and as the tree and it's extensions grows we can expand it to other pieces of software
making it accessible to all regardless of what you use. Alternatively, you could have a glance at the source code and copy/paste
it yourself.

```bash

pip install global-chem
pip install 'global-chem[graphing]'
pip install 'global-chem[forcefields]'
pip install 'global-chem[bioinformatics]'
pip install 'global-chem[cheminformatics]'
pip install 'global-chem[quantum_chemistry]'
pip install 'global-chem[development_operations]'
pip install 'global-chem[all]'

```

```python

from global_chem import GlobalChem
gc = GlobalChem()
gc.build_global_chem_network()
molecules = list(gc.get_node_smiles('pihkal').values())
print (molecules)
```

#### - Molecules Passing

|[![saythanks](https://img.shields.io/badge/RDKit-100%25-fg49b4.svg)](https://www.rdkit.org/)| [![saythanks](https://img.shields.io/badge/PartialSMILES-85.7%25-fg49b4.svg)](https://github.com/baoilleach/partialsmiles) | [![saythanks](https://img.shields.io/badge/DeepSMILES-99.25%25-lm89b5.svg)](https://github.com/baoilleach/deepsmiles) | [![saythanks](https://img.shields.io/badge/SELFIES-100%25-lm89b5.svg)](https://github.com/aspuru-guzik-group/selfies) | [![saythanks](https://img.shields.io/badge/MolVS-98.5%25-lm89b5.svg)](https://github.com/mcs07/MolVS) | [![saythanks](https://img.shields.io/badge/PySMILES-99.8%25-fg49b4.svg)](https://github.com/pckroon/pysmiles) |
|-|-|-|-|-|-|

GlobalChem Data Overview
=========================

### Nodes Contributors

Please follow the node contribution guidelines if you would like to elect your own or someone elses.

```
'global_chem': Node,                                                      # Suliman Sharif
'emerging_perfluoroalkyls': EmergingPerFluoroAlkyls,                      # Asuka Orr & Suliman Sharif
'montmorillonite_adsorption': MontmorilloniteAdsorption,                  # Asuka Orr & Suliman Sharif
'common_monomer_repeating_units': CommonMonomerRepeatingUnits,            # Suliman Sharif
'electrophilic_warheads_for_kinases': ElectrophilicWarheadsForKinases,    # Ruibin Liu & Suliman Sharif
'common_warheads_covalent_inhibitors': CommonWarheadsCovalentInhibitors,  # Shaoqi Zhan & Suliman Sharif
'rings_in_drugs': RingsInDrugs,                                           # Alexander Mackerell & Suliman Sharif
'iupac_blue_book_rings': IUPACBlueBookRings,                              # Suliman Sharif
'phase_2_hetereocyclic_rings': Phase2HetereoCyclicRings,                  # Suliman Sharif
'privileged_scaffolds': PrivilegedScaffolds,                              # Suliman Sharif
'iupac_blue_book': IUPACBlueBook,                                         # Suliman Sharif
'common_r_group_replacements': CommonRGroupReplacements,                  # Sunhwan Jo & Suliman Sharif
'braf_inhibitors': BRAFInhibitors,                                        # Aarion Romany & Suliman Sharif
'privileged_kinase_inhibitors': PrivilegedKinaseInhibitors,               # Suliman Sharif
'common_organic_solvents': CommonOrganicSolvents,                         # Suliman Sharif
'amino_acid_protecting_groups': AminoAcidProtectingGroups,                # Aziza Frank & Suliman Sharif
'schedule_one': ScheduleOne,                                              # Suliman Sharif
'schedule_two': ScheduleTwo,                                              # Suliman Sharif
'schedule_three': ScheduleThree,                                          # Suliman Sharif
'schedule_four': ScheduleFour,                                            # Suliman Sharif
'schedule_five': ScheduleFive,                                            # Suliman Sharif
'interstellar_space': InterstellarSpace,                                  # Suliman Sharif
'vitamins': Vitamins,                                                     # Suliman Sharif
'open_smiles': OpenSmiles,                                                # Suliman Sharif
'amino_acids': AminoAcids,                                                # Suliman Sharif
'pihkal': Pihkal,                                                         # Suliman Sharif
'nickel_ligands': NickelBidendatePhosphineLigands,                        # Suliman Sharif
'cimetidine_and_acyclovir': CimetidineAndAcyclovir,                       # Suliman Sharif
'common_regex_patterns': CommonRegexPatterns,                             # Chris Burke & Suliman Sharif
'how_to_live_longer': HowToLiveLonger,                                    # Suliman Sharif
'monoclonal_antibodies': MonoclonalAntibodies,                            # Asuka Orr & Suliman Sharif
'lube': Lube,                                                             # Daniel Khavrutskii & Suliman Sharif
'tainted_sexual_enhancements': TaintedSexualEnhancements,                 # Suliman Sharif
'exsens_products': ExsensProducts,                                        # Rebecca Pinette-Dorin & Suliman Sharif
'fda_list_one': FDAListOne,                                               # Mike Wostner & Suliman Sharif
'fda_list_two': FDAListTwo,                                               # Mike Wostner & Suliman Sharif
'fda_list_three': FDAListThree,                                           # Mike Wostner & Suliman Sharif
'fda_list_four': FDAListFour,                                             # Mike Wostner & Suliman Sharif
'fda_list_five': FDAListFive,                                             # Mike Wostner & Suliman Sharif
'fda_list_six': FDAListSix,                                               # Mike Wostner & Suliman Sharif
'fda_list_seven': FDAListSeven,                                           # Mike Wostner & Suliman Sharif
'constituents_of_cannabis_sativa': ConstituentsOfCannabisSativa,          # Ian Jones & Bettina Lier & Suliman Sharif
'phytocannabinoids': PhytoCannabinoids,                                   # Ian Jones & Bettina Lier & Suliman Sharif
'organophosphorous_nerve_agents': OrganoPhosphorousNerveAgents,           # Suliman Sharif
'organic_and_inorganic_bronsted_acids': OrganicAndInorganicBronstedAcids, # Nathaniel McClean & Suliman Sharif
'chemicals_from_biomass': ChemicalsFromBioMass,                           # Anthony Maiorana & Suliman Sharif
'salt': Salt,                                                             # Suliman Sharif
'drugs_from_snake_venom': DrugsFromSnakeVenom,                            # Suliman Sharif
'oral_contraceptives': OralContraceptives,                                # Suliman Sharif
'surfactants': Surfactants,                                               # Yiling Nan & Suliman Sharif
'lanthipeptides: LanthiPeptides                                           # Prabin Baral & Suliman Sharif
'alternative_jet_fuels': AlternativeJetFuels                              # Suliman Sharif
'mango_phytocompounds': MangoPhytoCompounds,                              # Damilola Bodun & Sevien Schulhoff & Suliman Sharif
'mango_amino_acids': MangoAminoAcids,                                     # Damilola Bodun & Sevien Schulhoff & Suliman Sharif
'mango_phenolic_acids': MangoPhenolicAcids,                               # Damilola Bodun & Sevien Schulhoff & Suliman Sharif
'mango_fatty_acids': MangoFattyAcids,                                     # Damilola Bodun & Sevien Schulhoff & Suliman Sharif
'mango_vitamins': MangoVitamins,                                          # Damilola Bodun & Sevien Schulhoff & Suliman Sharif
'mango_flavonoids': MangoFlavonoids,                                      # Damilola Bodun & Sevien Schulhoff & Suliman Sharif
'insect_sex_pheromones': InsectSexPheromones                              # Yuqing Liu & Suliman Sharif
```

| Chemical List                       | # of Entries | References                                                                                                                                                                                                                                                                                                           |
|-------------------------------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Amino Acids                         | 20           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Essential Vitamins                  | 13           | Common Knowledge                                                                                                                                                                                                                                                                                                     |
| Common Organic Solvents             | 42           | Fulmer, Gregory R., et al. “NMR Chemical Shifts of Trace Impurities: Common Laboratory Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist.”Organometallics, vol. 29, no. 9, May 2010, pp. 2176–79.                                                                          |
| Open Smiles                         | 94           | OpenSMILES Home Page. http://opensmiles.org/.                                                                                                                                                                                                                                                                        |
| IUPAC Blue Book (CRC Handbook) 2003 | 333          | Chemical Rubber Company. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edited by David R. Lide, 85. ed, CRC Press, 2004.                                                                                                                                               |
| Rings in Drugs                      | 92           | Taylor, Richard D., et al. “Rings in Drugs.” Journal of Medicinal Chemistry, vol. 57, no. 14, July 2014, pp. 5845–59. ACS Publications, https://doi.org/10.1021/jm4017625.                                                                                                                                           |
| Phase 2 Hetereocyclic Rings         | 19           | Broughton, Howard B., and Ian A. Watson. “Selection of Heterocycles for Drug Design.” Journal of Molecular Graphics & Modelling, vol. 23, no. 1, Sept. 2004, pp. 51–58. PubMed, https://doi.org/10.1016/j.jmgm.2004.03.016.                                                                                          |
| Privileged Scaffolds                | 47           | Welsch, Matthew E., et al. “Privileged Scaffolds for Library Design and Drug Discovery.” Current Opinion in Chemical Biology , vol. 14, no. 3, June 2010, pp. 347–61.PubMed, https://doi.org/10.1016/j.cbpa.2010.02.018.                                                                                             |
| Common Warheads                     | 29           | Gehringer, Matthias, and Stefan A. Laufer. “Emerging and Re-Emerging Warheads for Targeted Covalent Inhibitors: Applications in Medicinal Chemistry and Chemical Biology.”Journal of Medicinal Chemistry , vol. 62, no. 12, June 2019, pp. 5673–724. ACS Publications, https://doi.org/10.1021/acs.jmedchem.8b01153. |
| Common Polymer Repeating Units      | 78           | Hiorns, R. C., et al. “A brief guide to polymer nomenclature (IUPAC Technical Report).”Pure and Applied Chemistry , vol. 84, no. 10, Oct. 2012, pp. 2167–69., https://doi.org/10.1351/PAC-REP-12-03-05.                                                                                                              |
| Common R Group Replacements         | 499          | Takeuchi, Kosuke, et al. “R-Group Replacement Database for Medicinal Chemistry.”   Future Science OA , vol. 7, no. 8, Sept. 2021, p. FSO742.   future-science.com (Atypon) , https://doi.org/10.2144/fsoa-2021-0062.                                                                                                 |
| Electrophillic Warheads for Kinases | 24           | Petri, László, et al. “An Electrophilic Warhead Library for Mapping the Reactivity and Accessibility of Tractable Cysteines in Protein Kinases.” European Journal of Medicinal Chemistry, vol. 207, Dec. 2020, p. 112836. PubMed, https://doi.org/10.1016/j.ejmech.2020.112836.                                      |
| Privileged Scaffolds for Kinases    | 29           | Hu, Huabin, et al. “Systematic Comparison of Competitive and Allosteric Kinase Inhibitors Reveals Common Structural Characteristics.” European Journal of Medicinal Chemistry, vol. 214, Mar. 2021, p. 113206. ScienceDirect, https://doi.org/10.1016/j.ejmech.2021.113206.                                          |
| BRaf Inhibitors                     | 54           | Agianian, Bogos, and Evripidis Gavathiotis. “Current Insights of BRAF Inhibitors in Cancer.” Journal of Medicinal Chemistry, vol. 61, no. 14, July 2018, pp. 5775–93. ACS Publications, https://doi.org/10.1021/acs.jmedchem.7b01306.                                                                                |
| Common Amino Acid Protecting Groups | 346          | Isidro-Llobet, Albert, et al. “Amino Acid-Protecting Groups.” Chemical Reviews, vol. 109, no. 6, June 2009, pp. 2455–504. DOI.org (Crossref), https://doi.org/10.1021/cr800323s.                                                                                                                                     |
| Emerging Perfluoroalkyls            | 27           | Pelch, Katherine E., et al. “PFAS Health Effects Database: Protocol for a Systematic Evidence Map.” Environment International, vol. 130, Sept. 2019, p. 104851. ScienceDirect, https://doi.org/10.1016/j.envint.2019.05.045.                                                                                         |
| Chemicals For Clay Adsorption       | 33           | Orr, Asuka A., et al. “Combining Experimental Isotherms, Minimalistic Simulations, and a Model to Understand and Predict Chemical Adsorption onto Montmorillonite Clays.” ACS Omega, vol. 6, no. 22, June 2021, pp. 14090–103. PubMed, https://doi.org/10.1021/acsomega.1c00481.                                     |
| Schedule 1 United States Narcotics  | 240          | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 2 United States Narcotics  | 60           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 3 United States Narcotics  | 22           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 4 United States Narcotics  | 77           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 5 United States Narcotics  | 8            | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Pihkal                              | 179          | Shulgin, Alexander T., and Ann Shulgin. Pihkal: A Chemical Love Story. 1. ed., 8. print, Transform, 2010.                                                                                                                                                                                                            |
| Excipients Cimetidine & Acyclovir   | 14           | Vaithianathan, Soundarya, et al. “Effect of Common Excipients on the Oral Drug Absorption of Biopharmaceutics Classification System Class 3 Drugs Cimetidine and Acyclovir.” Journal of Pharmaceutical Sciences, vol. 105, no. 2, Feb. 2016, pp. 996–1005. PubMed, https://doi.org/10.1002/jps.24643.                |
| Nickel Bidendate Phosphine Ligands  | N/A          | Clevenger, Andrew L., et al. “Trends in the Usage of Bidentate Phosphines as Ligands in Nickel Catalysis.” Chemical Reviews, vol. 120, no. 13, July 2020, pp. 6124–96. DOI.org (Crossref), https://doi.org/10.1021/acs.chemrev.9b00682.                                                                              |
| HowToLiveLonger                     | 4            | https://github.com/geekan/HowToLiveLonger                                                                                                                                                                                                                                                                            |
| Monoclonal Antibodies               | 19           | https://labels.fda.gov/                                                                                                                                                                                                                                                                                              |
| Common Lubricants for Sex           | 38           | https://exsens-usa.com/blogs/your-body-your-pleasure/lube-lessons-glossary-of-common-sex-lube-ingredients                                                                                                                                                                                                            |
| Tainted Sexual Enhancements         | 4            | FDA Tainted Sexual Enhancements                                                                                                                                                                                                                                                                                      |
| Salt                                | 14           | OpenFoodFacts https://github.com/openfoodfacts                                                                                                                                                                                                                                                                       |
| Exsens Sexual Wellness              | 59           | https://exsens-usa.com/                                                                                                                                                                                                                                                                                              |
| FDA Color Additive List 1           | 12           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| FDA Color Additive List 2           | 15           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| FDA Color Additive List 3           | 16           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| FDA Color Additive List 4           | 39           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| FDA Color Additive List 5           | 27           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| FDA Color Additive List 6           | 29           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| FDA Color Additive List 7           | 37           | https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list                                                                                                                                                                                                                                   |
| Constituents of Cannabis Sativa     | 394          | Turner, C. E., et al. “Constituents of Cannabis Sativa L. XVII. A Review of the Natural Constituents.” Journal of Natural Products, vol. 43, no. 2, Apr. 1980, pp. 169–234. PubMed                                                                                                                                   |
| Phytocannabinoids                   | 111          | Hanuš, Lumír Ondřej, et al. “Phytocannabinoids: A Unified Critical Inventory.” Natural Product Reports, vol. 33, no. 12, Nov. 2016, pp. 1357–92. PubMed,                                                                                                                                                             |
| OrganoPhosphorous Nerve Agents      | 14           | Mukherjee, Sudisha, and Rinkoo Devi Gupta. “Organophosphorus Nerve Agents: Types, Toxicity, and Treatments.” Journal of Toxicology, vol. 2020, Sept. 2020, p. 3007984.                                                                                                                                               |
| Cengage Bronsted Acids              | 42           | https://cxp.cengage.com/contentservice/assets/owms01h/references/chemtables/org_chem/pKaTable.html                                                                                                                                                                                                                   |
| Chemicals From Biomass              | 17           | Wittcoff, Harold A., et al. Industrial Organic Chemicals: Wittcoff/Organic Chemicals. John Wiley & Sons, Inc., 2004                                                                                                                                                                                                  |
| Drugs From Snake Venom              | 7            | Oliveira, Ana L., et al. “The Chemistry of Snake Venom and Its Medicinal Potential.” Nature Reviews Chemistry, vol. 6, no. 7, July 2022, pp. 451–69                                                                                                                                                                  |
| Oral Contraceptives                 | 17           | Coleman, William F. “The Molecules of Oral Contraceptives.” Journal of Chemical Education, vol. 87, no. 7, July 2010, pp. 760–61.                                                                                                                                                                                    |
| Surfactants for Skin                | 36           | Date, Abhijit A., and Vandana B. Patravale. “Microemulsions: Applications in Transdermal and Dermal Delivery.” Critical Reviews&trade; in Therapeutic Drug Carrier Systems, vol. 24, no. 6, 2007.                                                                                                                    |
| LanthiPeptides                      | 2            | Pokhrel, Rudramani, et al. “Molecular Mechanisms of Pore Formation and Membrane Disruption by the Antimicrobial Lantibiotic Peptide Mutacin 1140.” Physical Chemistry Chemical Physics, vol. 21, no. 23, June 2019, pp. 12530–39.                                                                                    |
| Alternative Jet Fuels               | 59           | Chemical Composition and Fuel Properties of Alternative Jet Fuels :: BioResources. https://bioresources.cnr.ncsu.edu/.                                                                                                                                                                                               |
| Mango Amino Acids                   | 19           | Maldonado-Celis, Maria Elena, et al. “Chemical Composition of Mango (Mangifera Indica L.) Fruit: Nutritional and Phytochemical Compounds.” Frontiers in Plant Science, vol. 10, Oct. 2019, p. 1073.                                                                                                                  |
| Mango Phenoloic Acids               | 10           | Maldonado-Celis, Maria Elena, et al. “Chemical Composition of Mango (Mangifera Indica L.) Fruit: Nutritional and Phytochemical Compounds.” Frontiers in Plant Science, vol. 10, Oct. 2019, p. 1073.                                                                                                                  |
| Mango Fatty Acids                   | 24           | Maldonado-Celis, Maria Elena, et al. “Chemical Composition of Mango (Mangifera Indica L.) Fruit: Nutritional and Phytochemical Compounds.” Frontiers in Plant Science, vol. 10, Oct. 2019, p. 1073.                                                                                                                  |
| Mango Vitamins                      | 10           | Maldonado-Celis, Maria Elena, et al. “Chemical Composition of Mango (Mangifera Indica L.) Fruit: Nutritional and Phytochemical Compounds.” Frontiers in Plant Science, vol. 10, Oct. 2019, p. 1073.                                                                                                                  |
| Mango Flavonoids                    | 11           | Maldonado-Celis, Maria Elena, et al. “Chemical Composition of Mango (Mangifera Indica L.) Fruit: Nutritional and Phytochemical Compounds.” Frontiers in Plant Science, vol. 10, Oct. 2019, p. 1073.                                                                                                                  |
| Insect Sex Pheromones               | 37           | Jacobson, Martin. Insect Sex Pheromones. New York, Academic Press, 1992.                                                      |

Features
========

| Extension                       | Description                                                                                                             | Appplication   |
|---------------------------------|-------------------------------------------------------------------------------------------------------------------------|------------------|
| GlobalChem Chemical Entities    | GlobalChem has internal Molecule objects with all common attributes associated and conversion to SMILES                 | forcefields       |
| GlobalChem Biological Entities  | GlobalChem has internal DNA/RNA/Protein/Molecule objects with all common attributes associated and conversion to SMILES | bioinformatics   |
| Visualize DNA/RNA Strands       | Visualize DNA and RNA Strands and add labels to them | bioinformatics   |
| ForceField Molecules            | GlobalChem can parse, manipulate, and write CGenFF and GaFF2 files as objects                                           | forcefields       |
| PDF Generation and Parsing      | GlobalChem can generate SMILES to PDF and convert the PDF to SMILES                                                     | cheminformatics              |
| SMILES Validation               | GlobalChem has connection to PySMILES, DeepSMILES, PartialSmiles, SELFIES, MolVS for validation of SMILES sets          | cheminformatics       |
| SMILES Protonation States       | GlobalChem can take a set of compounds and predict the protonation states of a SMILES string over a range of pH         | chemfinformatics       |
| Open Source Database Monitoring | GlobalChem uses Uptime-Cheminformatics to Keep Track of Open Source Chemical Data                                       | development_operations       |
| Networkx Software Adapter       | GlobalChem Network can be converted into NetworkX Graph Objects                                                         | cheminformatics       |
| SMARTS Pattern Validation       | GlobalChem uses the MiniFrag Database to test SMARTS strings accuracy for functional group selection                    | cheminformatics       |
| Principal Component Analysis    | GlobalChem can readily interpret SMILES, fingerprint, cluster and apply PCA analysis user can tweak parameters          | cheminformatics |
| Drug Design Filters             | GlobalChem can filter compounds based on Common Drug Design Filtering Rules                                            | cheminformatics       |
| Deep Layer Scatter Analysis     | To visualize relations between sets of molecules, GlobalChem offers a parallel coordinate diagram generation            | cheminformatics | 
| Sunbursting Radial Analysis     | GlobalChem offers a sunbursting mechanism to allow uses to observe how sets of compounds relate to the common set      | cheminformatics |
| Graphing Templates              | GlobalChem offers graphing templates to aid in faster data analysis, currently the only offer is Plotly               | cheminformatics |
| CGenFF Dissimilarity Score      | GlobalChem can offer the difference between two molecules based on their Atom Types                                     | forcefields       |
| OneHot Encoding                 | GlobalChem has it's own one hot encoder and decoder based on the common lists for Machine Learning                      | cheminformatics       |
| SMARTS Pattern Identifier       | GlobalChem connects to the SMARTS Plus and can offer visualization into different SMARTS components                     | cheminformatics       |
| Psi4 Parser       | Offer parsing of Psi4 Output Files and extracting values                    | quantum_chemistry  |
| Coordinate Store       | A warehouse for coodinates of small molecules for distribution in xyz and zm-matrix                   | quantum_chemistry  |
| Visualize Molecular Orbitals       | Visualize the Cube Files from Psi4 Output cubeprop                   | quantum_chemistry  |
‌
Contributions
=============

![Alt](https://repobeats.axiom.co/api/embed/a5562010add209dab03a109cd9feaef6234dd8ca.svg "Repobeats analytics image")

Licensing
=========

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)
