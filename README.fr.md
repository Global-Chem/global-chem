<h1 align="center">Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities </h1>

Global Chem est une base de données de graphes open source et une API pour les listes de produits chimiques courants et rares utilisant IUPAC en entrée et SMILES/SMARTS en sortie. Comme
dont j'ai surtout besoin pour moi-même alors que je cherche à travers l'infini chimique.

J'ai trouvé ces listes écrites dans l'histoire utiles, elles proviennent de différents domaines mais sont agrégées
dans le format le plus courant des chimistes organiques (IUPAC) et le langage commun du cheminformaticien (SMILES) et pour
correspondance de motifs (SMART).

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

# Aperçu

## GlobalChem

#### - Commencer

| [![Documentation](https://img.shields.io/badge/GitBook-Docu-lightblue)](https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/) | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1HaDvAYpaX0_2Eqx1NSl5uJTRkE_0na-L?usp=sharing) | [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.chemicalgraphtheory.com) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) |
| -------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |

#### - Validation

| [![saythanks](https://img.shields.io/badge/RDKit-100%25-fg49b4.svg)](https://www.rdkit.org/) | [![saythanks](https://img.shields.io/badge/PartialSMILES-85.7%25-fg49b4.svg)](https://github.com/baoilleach/partialsmiles) | [![saythanks](https://img.shields.io/badge/DeepSMILES-99.25%25-lm89b5.svg)](https://github.com/baoilleach/deepsmiles) | [![saythanks](https://img.shields.io/badge/SELFIES-100%25-lm89b5.svg)](https://github.com/aspuru-guzik-group/selfies) | [![saythanks](https://img.shields.io/badge/MolVS-98.5%25-lm89b5.svg)](https://github.com/mcs07/MolVS) | [![saythanks](https://img.shields.io/badge/PySMILES-99.8%25-fg49b4.svg)](https://github.com/pckroon/pysmiles) |
| -------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |

#### - Statistiques

| [![Downloads](https://pepy.tech/badge/global-chem)](https://pepy.tech/project/global-chem) | ![visitor badge](https://visitor-badge.glitch.me/badge?page_id=sulstice.global-chem) | [![Man Hours](https://img.shields.io/endpoint?url=https%3A%2F%2Fmh.jessemillar.com%2Fhours%3Frepo%3Dhttps%3A%2F%2Fgithub.com%2FSulstice%2Fglobal-chem)](https://jessemillar.com/r/man-hours) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) |
| ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |

#### - Actions Github

| [![Test System](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml) | [![pre-commit.ci status](https://results.pre-commit.ci/badge/github/Sulstice/global-chem/master.svg)](https://results.pre-commit.ci/latest/github/Sulstice/global-chem/master) | [![publish](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml) | [![Translate README](https://github.com/Sulstice/global-chem/actions/workflows/translate_readme.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/translate_readme.yml) |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

#### - Informations sur la construction

| [![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem) | [![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) | [![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250) | [![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_shield) | [![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/) | [![Maturity badge - level 2](https://img.shields.io/badge/Maturity-Level%202%20--%20First%20Release-yellowgreen.svg)](https://github.com/tophat/getting-started/blob/master/scorecard.md) |
| ------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

#### - GlobalChemExtensions

| [![PyPI version](https://badge.fury.io/py/global-chem-extensions.svg)](https://badge.fury.io/py/global-chem-extensions) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) | [![Downloads](https://pepy.tech/badge/global-chem-extensions)](https://pepy.tech/project/global-chem-extensions) |
| ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- |

### Bande démo d'extension

| Application                 | Introduction                                                                                                                                               | Utilisation avancée                                                                                 |                                                                                                                                                            |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| champs de force             | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1hW0K6V5zPDFdZvFkYarbOr9wRoof2n4s?usp=sharing) | CGenFF Molecule Loader et similarité des types d'atomes                                             | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1OEm1dACd_wkw_JQITemy18pgfKdE4XlX?usp=sharing) |
| bioinformatique             | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sCNw9FIcFfTEgmrdokADXgJGnsPEOGbX?usp=sharing) | Utilisez l'algorithme de Bostrom pour filtrer les ligands par PDB                                   | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1a3Ys1rpqFzxBz95DQwJVHfunBxpywP7f?usp=sharing) |
| cheminformatics             | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1z0ilrakoRJ8maapMNHwtPf83pKK1THST?usp=sharing) |                                                                                                     |                                                                                                                                                            |
| chimie quantique            | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1BGLQphP1IMLndeyavHM_6_qXQJI_7gFU?usp=sharing) | Tracer la théorie quantique et l'ensemble de base par rapport à l'hamiltonien des petites molécules | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/11JGsen912TmyR5Dds-Opon5cRkrtVoT7?usp=sharing) |
| opérations_de_développement | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1v_yXPXPilbWUZGkel_yekBr3EFnIUNUL?usp=sharing) |                                                                                                     |                                                                                                                                                            |

# Installation

GlobalChem va être distribué via PyPi et au fur et à mesure que l'arborescence et ses extensions grandissent, nous pouvons l'étendre à d'autres logiciels
le rendant accessible à tous, peu importe ce que vous utilisez. Alternativement, vous pouvez jeter un coup d'œil au code source et copier/coller
le toi-même.

```bash

pip install global-chem

```

Si vous souhaitez installer le package d'extensions pour des fonctionnalités supplémentaires, chaque application peut être installée indépendamment l'une de l'autre ou vous pouvez toutes les installer avec le`all`:

```bash

pip install 'global-chem[graphing]'
pip install 'global-chem[forcefields]'
pip install 'global-chem[bioinformatics]'
pip install 'global-chem[cheminformaitcs]'
pip install 'global-chem[quantum_chemistry]'
pip install 'global-chem[development_operations]'
pip install 'global-chem[all]'

```

# Démarrage rapide

Ici, nous chargeons le`global-chem[cheminformatics]`paquet d'extensions et le`GlobalChem`arbre. Nous extrayons SMILES du livre populaire, Pihkal, et effectuons une analyse en composantes principales chimio-informatique sur la liste chimique.

```python

from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
gc_cheminfo = GlobalChemExtensions().cheminformatics()

gc.build_global_chem_network()
smiles_list = list(gc.get_node_smiles('pihkal').values())

print (f"SMILES: {smiles_list[0]}")

gc_cheminfo.node_pca_analysis(smiles_list)

```

# GlobalChem

### Règles

Le Graph Network (GN) est livré avec quelques règles qui, pour l'instant, facilitent l'ingénierie logicielle pour le développeur.

-   Il doit y avoir un nœud racine.
-   Lors de l'ajout d'un nœud, chaque nœud doit être connecté.
-   Pour supprimer un nœud, il ne doit pas avoir d'enfant.

Le Deep Graph Network (DGN) est également livré avec quelques règles pour faciliter la mise en œuvre :

-   Il doit y avoir un nœud racine de 1 qui marque comme votre nœud "d'entrée".
-   Lors de l'ajout d'une couche, tous les nœuds seront ajoutés à toutes les couches précédentes en tant qu'enfants. (Les gens peuvent utiliser la fonction de suppression de nœud pour effectuer des abandons).

### Graphique des connaissances

Juste sans dépendances, initialisez la classe et c'est parti ! Tous les groupes communs et rares du monde
a ta disposition.

<p align="center">
  <img width="850" height="650" src="images/figures/figure_1_new.png">
</p>

### Nœuds Contributeurs

Veuillez suivre les directives de contribution des nœuds si vous souhaitez élire vous-même ou quelqu'un d'autre.

    'global_chem': Node,
    'emerging_perfluoroalkyls': EmergingPerFluoroAlkyls,                      # Asuka Orr & Suliman Sharif
    'montmorillonite_adsorption': MontmorilloniteAdsorption,                  # Asuka Orr & Suliman Sharif
    'common_monomer_repeating_units': CommonMonomerRepeatingUnits,            # Suliman Sharif
    'electrophilic_warheads_for_kinases': ElectrophilicWarheadsForKinases,    # Ruibin Liu & Suliman Sharif
    'common_warheads_covalent_inhibitors': CommonWarheadsCovalentInhibitors,  # Shaoqi Zhao & Suliman Sharif
    'rings_in_drugs': RingsInDrugs,                                           # Alexander Mackerell & Suliman Sharif
    'iupac_blue_book_rings': IUPACBlueBookRings,                              # Suliman Sharif
    'phase_2_hetereocyclic_rings': Phase2HetereoCyclicRings,                  # Suliman Sharif
    'privileged_scaffolds': PrivilegedScaffolds,                              # Suliman Sharif
    'iupac_blue_book': IUPACBlueBook,                                         # Suliman Sharif
    'common_r_group_replacements': CommonRGroupReplacements,                  # Sunhwan Jo & Suliman Sharif
    'braf_inhibitors': BRAFInhibitors,                                        # Aarion Romany
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

<details><summary><h3>Node List</h1><br/></summary>

| Liste des produits chimiques                  | # d'entrées | Références                                                                                                                                                                                                                                                                                                                                |
| --------------------------------------------- | ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Acides aminés                                 | 20          | Connaissance commune                                                                                                                                                                                                                                                                                                                      |
| Vitamines essentielles                        | 13          | Connaissance commune                                                                                                                                                                                                                                                                                                                      |
| Solvants organiques courants                  | 42          | Fulmer, Gregory R., et al. "Déplacements chimiques RMN des traces d'impuretés : solvants de laboratoire courants, composés organiques et gaz dans les solvants deutérés pertinents pour le chimiste organométallique." Organométalliques, vol. 29, non. 9, mai 2010, p. 2176–79.                                                          |
| Sourires ouverts                              | 94          | Page d'accueil d'OpenSMILES.<http://opensmiles.org/>.                                                                                                                                                                                                                                                                                     |
| Livre bleu IUPAC (manuel CRC) 2003            | 333         | Société de caoutchouc chimique. CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edité par David R. Lide, 85. ed, CRC Press, 2004.                                                                                                                                                             |
| Anneaux dans la drogue                        | 92          | Taylor, Richard D., et al. "Anneaux dans la drogue." Tourillon de chimie médicinale, vol. 57, non. 14, juillet 2014, p. 5845–59. AEC Publications,<https://doi.org/10.1021/jm4017625>.                                                                                                                                                    |
| Anneaux hétérocycliques de phase 2            | 19          | Broughton, Howard B., et Ian A. Watson. "Sélection d'hétérocycles pour la conception de médicaments." Journal of Molecular Graphics & Modeling, vol. 23, non. 1, septembre 2004, p. 51–58. PubMed,<https://doi.org/10.1016/j.jmgm.2004.03.016>.                                                                                           |
| Échafaudages privilégiés                      | 47          | Welsch, Matthew E., et al. "Échafaudages privilégiés pour la conception de bibliothèques et la découverte de médicaments." Opinion actuelle en biologie chimique, vol. 14, non. 3, juin 2010, p. 347–61. PubMed,<https://doi.org/10.1016/j.cbpa.2010.02.018>.                                                                             |
| Ogives communes                               | 29          | Gehringer, Matthias et Stefan A. Laufer. "Ogives émergentes et réémergentes pour les inhibiteurs covalents ciblés : applications en chimie médicinale et en biologie chimique."Journal of Medicinal Chemistry, vol. 62, non. 12, juin 2019, p. 5673–724. AEC Publications,<https://doi.org/10.1021/acs.jmedchem.8b01153>.                 |
| Unités répétitives polymères communes         | 78          | Hiorns, R.C., et al. "Un bref guide de la nomenclature des polymères (rapport technique IUPAC)." Pure and Applied Chemistry, vol. 84, non. 10, oct. 2012, p. 2167–69.,<https://doi.org/10.1351/PAC-REP-12-03-05>.                                                                                                                         |
| Remplacements courants du groupe R            | 499         | Takeuchi, Kosuke et al. "Base de données de remplacement du groupe R pour la chimie médicinale." Future Science OA, vol. 7, non. 8, septembre 2021, p. OFS742. future-science.com (Atypon) ,<https://doi.org/10.2144/fsoa-2021-0062>.                                                                                                     |
| Têtes électrophiles pour kinases              | 24          | Petri, Laszlo, et al. "Une bibliothèque d'ogives électrophiles pour cartographier la réactivité et l'accessibilité des cystéines traitables dans les protéines kinases." Journal européen de chimie médicinale, vol. 207, décembre 2020, p. 112836. PubMed,<https://doi.org/10.1016/j.ejmech.2020.112836>.                                |
| Échafaudages privilégiés pour les kinases     | 29          | Hu, Huabin et al. "La comparaison systématique des inhibiteurs de kinase compétitifs et allostériques révèle des caractéristiques structurelles communes." Journal européen de chimie médicinale, vol. 214, mars 2021, p. 113206. ScienceDirect,<https://doi.org/10.1016/j.ejmech.2021.113206>.                                           |
| Inhibiteurs de BRAF                           | 54          | Agianian, Bogos et Evripidis Gavathiotis. « Aperçus actuels des inhibiteurs de BRAF dans le cancer. Tourillon de chimie médicinale, vol. 61, non. 14, juillet 2018, p. 5775–93. AEC Publications,<https://doi.org/10.1021/acs.jmedchem.7b01306>.                                                                                          |
| Groupes protecteurs d'acides aminés courants  | 346         | Isidro-Llobet, Albert, et al. "Groupes de protection des acides aminés." Revues chimiques, vol. 109, non. 6, juin 2009, p. 2455–504. DOI.org (référence croisée),<https://doi.org/10.1021/cr800323s>.                                                                                                                                     |
| Perfluoroalkyles émergents                    | 27          | Pelch, Katherine E., et al. "Base de données des effets sur la santé des PFAS : protocole pour une carte de preuves systématiques." Environnement International, vol. 130, septembre 2019, p. 104851. ScienceDirect,<https://doi.org/10.1016/j.envint.2019.05.045>.                                                                       |
| Produits chimiques pour l'adsorption d'argile | 33          | Orr, Asuka A., et al. "Combinant des isothermes expérimentales, des simulations minimalistes et un modèle pour comprendre et prédire l'adsorption chimique sur les argiles de montmorillonite." ACS Oméga, vol. 6, non. 22, juin 2021, p. 14090–103. PubMed,<https://doi.org/10.1021/acsomega.1c00481>.                                   |
| Cannabinoïdes                                 | 63          | Turner, Carlton E., et al. “Constituants du Cannabis Sativa L. XVII. Un examen des constituants naturels. Journal des produits naturels, vol. 43, non. 2, mars 1980, p. 169–234. AEC Publications,<https://doi.org/10.1021/np50008a001>.                                                                                                  |
| Annexe 1 États-Unis Stupéfiants               | 240         | ECFR :: 21 CFR Partie 1308 - Annexes.                                                                                                                                                                                                                                                                                                     |
| Annexe 2 États-Unis Stupéfiants               | 60          | ECFR :: 21 CFR Partie 1308 - Annexes.                                                                                                                                                                                                                                                                                                     |
| Annexe 3 États-Unis Stupéfiants               | 22          | ECFR :: 21 CFR Partie 1308 - Annexes.                                                                                                                                                                                                                                                                                                     |
| Annexe 4 États-Unis Stupéfiants               | 77          | ECFR :: 21 CFR Partie 1308 - Annexes.                                                                                                                                                                                                                                                                                                     |
| Annexe 5 États-Unis Stupéfiants               | 8           | ECFR :: 21 CFR Partie 1308 - Annexes.                                                                                                                                                                                                                                                                                                     |
| Pihkal                                        | 179         | Shulgin, Alexander T., et Ann Shulgin. Pihkal : Une histoire d'amour chimique. 1. éd., 8. impression, Transformer, 2010.                                                                                                                                                                                                                  |
| Excipients Cimétidine & Acyclovir             | 14          | Vaithianathan, Soundarya, et al. "Effet des excipients communs sur l'absorption orale des médicaments du système de classification biopharmaceutique des médicaments de classe 3 cimétidine et acyclovir." Journal des sciences pharmaceutiques, vol. 105, non. 2, février 2016, p. 996–1005. PubMed,<https://doi.org/10.1002/jps.24643>. |
| Ligands de phosphine bidendate de nickel      | N/A         | Clevenger, Andrew L., et al. "Tendances dans l'utilisation des phosphines bidentées comme ligands dans la catalyse au nickel." Revues chimiques, vol. 120, non. 13, juillet 2020, p. 6124–96. DOI.org (référence croisée),<https://doi.org/10.1021/acs.chemrev.9b00682>.                                                                  |
| Modèles de regex courants                     | 1           |                                                                                                                                                                                                                                                                                                                                           |

</details>                            

# GlobalChemExtensions

Une variété d'outils sont disponibles pour vous permettre de parcourir et d'analyser les données et la liste complète des différentes applications peut être trouvée dans la démo google colab ou la documentation Gitbook. Une démonstration des extensions de visualisation de données conçues avec plotly et bokeh est présentée ci-dessous :

<p align="center">
  <img width="800" height="600" src="https://raw.githubusercontent.com/Sulstice/global-chem/master/images/figures/figure_10.png">
</p>

<details><summary><h3>Extension List</h1><br/></summary>

| Extension                                      | La description                                                                                                                                                   | Application                 |
| ---------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| Entités chimiques GlobalChem                   | GlobalChem a des objets Molecule internes avec tous les attributs communs associés et la conversion en SMILES                                                    | champs de force             |
| Entités biologiques GlobalChem                 | GlobalChem a des objets internes d'ADN/ARN/protéine/molécule avec tous les attributs communs associés et la conversion en SOURIRES                               | bioinformatique             |
| Visualiser les brins d'ADN/ARN                 | Visualisez les brins d'ADN et d'ARN et ajoutez-leur des étiquettes                                                                                               | bioinformatique             |
| Molécules de champ de force                    | GlobalChem peut analyser, manipuler et écrire des fichiers CGenFF et GaFF2 en tant qu'objets                                                                     | champs de force             |
| Génération et analyse de PDF                   | GlobalChem peut générer des SMILES en PDF et convertir le PDF en SMILES                                                                                          | cheminformatics             |
| Validation SOURIRES                            | GlobalChem est connecté à PySMILES, DeepSMILES, PartialSmiles, SELFIES, MolVS pour la validation des ensembles SMILES                                            | cheminformatics             |
| SOURIRES États de protonation                  | GlobalChem peut prendre un ensemble de composés et prédire les états de protonation d'une chaîne SMILES sur une plage de pH                                      | chemininformatique          |
| Surveillance de la base de données open source | GlobalChem utilise Uptime-Cheminformatics pour suivre les données chimiques open source                                                                          | opérations_de_développement |
| Adaptateur logiciel Networkx                   | Le réseau GlobalChem peut être converti en objets graphiques NetworkX                                                                                            | cheminformatics             |
| Validation du modèle SMARTS                    | GlobalChem utilise la base de données MiniFrag pour tester la précision des chaînes SMARTS pour la sélection des groupes fonctionnels                            | cheminformatics             |
| Analyse des composants principaux              | GlobalChem peut facilement interpréter les SMILES, les empreintes digitales, les clusters et appliquer l'analyse PCA, l'utilisateur peut modifier les paramètres | cheminformatics             |
| Filtres de conception de médicaments           | GlobalChem peut filtrer les composés en fonction des règles communes de filtrage de la conception des médicaments                                                | cheminformatics             |
| Analyse de diffusion en couche profonde        | Pour visualiser les relations entre des ensembles de molécules, GlobalChem propose une génération de diagrammes de coordonnées parallèles                        | cheminformatics             |
| Analyse radiale Sunbursting                    | GlobalChem propose un mécanisme de sunbursting pour permettre aux utilisateurs d'observer comment des ensembles de composés sont liés à l'ensemble commun        | cheminformatics             |
| Modèles graphiques                             | GlobalChem propose des modèles de graphiques pour faciliter une analyse plus rapide des données, actuellement la seule offre est Plotly                          | cheminformatics             |
| Score de dissimilitude CGenFF                  | GlobalChem peut offrir la différence entre deux molécules en fonction de leurs types d'atomes                                                                    | champs de force             |
| Un encodage à chaud                            | GlobalChem possède son propre encodeur et décodeur à chaud basé sur les listes communes pour l'apprentissage automatique                                         | cheminformatics             |
| Identificateur de modèle SMARTS                | GlobalChem se connecte au SMARTS Plus et peut offrir une visualisation dans différents composants SMARTS                                                         | cheminformatics             |
| Analyseur Psi4                                 | Offre l'analyse des fichiers de sortie Psi4 et l'extraction des valeurs                                                                                          | chimie quantique            |
| Magasin de coordonnées                         | Un entrepôt de coordonnées de petites molécules pour distribution en xyz et z-matrix                                                                             | chimie quantique            |
| Visualiser les orbitales moléculaires          | Visualisez les fichiers de cube à partir de Psi4 Output cubeprop                                                                                                 | chimie quantique            |

</details>

# Conformité des logiciels open source

`GlobalChem`suit les mêmes principes énoncés dans la partie 11 du titre 21 du Code of Federal Regulations ; Dossiers électroniques,
Document d'orientation sur les signatures électroniques (21 CFR Part 11). Puisqu'il n'y a pas de directives formelles sur la façon dont les logiciels open source doivent être manipulés, nous
tenter de répondre aux exigences. La FDA considère que la partie 11 s'applique aux critères suivants d'enregistrements électroniques et comment`GlobalChem`accomplit chaque composant :

-   **Plausibilité :**`GlobalChem`a été construit sur des données extraites de livres et d'articles en utilisant la lecture et le redessin. Il ajoute une composante de
    Chaînes IUPAC/SMILES/SMARTS pour le stocker électroniquement, ce qui donne à ses données son composant unique. Les enregistrements sont open source
    et la version contrôlée de manière appropriée par les mainteneurs du référentiel et les commentaires de la communauté open source.`GlobalChem`Les objectifs de sont encore inconnus car il entre dans le déploiement open source. Nous avons construit des fonctions étendues qui vivent dans
    un paquet séparé`GlobalChemExtensions`qui dépendent de`GlobalChem`. Étant donné que chaque version est emballée de manière appropriée, de
    s'appuyer sur une version est un besoin alors son logiciel est disponible sur`Github`et`PyPi`. Une Procédure Opérationnelle Standard (SOP)
    peut être déposé à partir de la documentation de l'utilitaire d'extensions conservée sur`Gitbook`.

-   **Validation:**`GlobalChem`respecte la catégorie 3 des bonnes pratiques de fabrication automatisées (GAMP) qui correspond à "un logiciel utilisé tel qu'il est installé"
    et potentiellement "configurable".`GlobalChem`les tests viennent de l'intérieur, la documentation sert de test ultime
    pour la fonctionnalité car c'est ce que les utilisateurs testeront le plus puisque nous nous appuyons sur l'open source. Une intégration continue (IC)
    système est également construit simultanément pour servir de test de fonctionnalité de base du`GlobalChem`réseau de graphes. Les Données stockées
    est maintenu par des experts dans le domaine mais peut être modifié en fonction des commentaires de la communauté si une erreur est détectée.

-   **Piste d'audit :**`GlobalChem`est sous contrôle de version avec`Git`et hébergé sur la plateforme de Microsoft`Github`.`GlobalChem`suit une sémantique
    contrôle de version du schéma`X1.X2.X3`:`X1`marque les versions stables formelles avec des tests et une documentation et signifie
    gros refactoring au logiciel ou en fonctionnalité,`X2`signifie qu'une nouvelle fonctionnalité est ajoutée avec ou sans tests et documentation, mais
    itère comme tel.`X3`signifie un correctif "à chaud" (quelque chose qui est un bogue facile), une petite fonctionnalité ou un paramètre supplémentaire à ajouter à une fonction
    , ou itération aux données.

-   **Systèmes hérités :**`GlobalChem`est opérationnel depuis près de 2 ans depuis sa première version avec la version`0.3.0`en mai 2020.`GlobalChem`a été construit avec un parcours complet dans la communauté open source avec chaque version cataloguée et une visibilité pour tous. Cela satisfait
    les grandes lignes des règles pour déterminer un système hérité. Nous utilisons les commentaires de la communauté fournis par les plateformes de médias sociaux (Twitter, Github, LinkedIn)
    comme preuve documentée et justification que`GlobalChem`est adapté à l'utilisation prévue de cheminformatics.

-   **Copies des enregistrements :**`GlobalChem`a des enregistrements stockés sur`Github`pour le logiciel qui peut être exporté vers une variété de formats fournis par
    Microsoft. Pour la documentation, il est hébergé sur`Gitbook`et version contrôlée conformément au logiciel. Chaque "livre"
    peut être exporté au format de données portable (PDF) approprié pour la soumission à la FDA.

-   **Conservation des dossiers:**`GlobalChem`a un enregistrement de la documentation versionnée contrôlée par un identifiant unique (UUID) qui sert d'identifiant
    pour chaque itération stockée sur`Gitbook`. Chaque version est stockée sous forme de fichiers Markdown et peut être convertie en PDF, si nécessaire.

`GlobalChem`possède une licence publique Mozilla version 2.0.`GlobalChem`vous permet d'utiliser le logiciel dans votre travail plus large et
étendez-le avec des modifications si vous le souhaitez. L'éventualité est que si vous installez`GlobalChem`et publier un nouveau logiciel
alors vous devez suivre les mêmes principes installés dans notre licence pour la communauté open source.

# Collecte de données

Les références et les listes de composés associés sont sélectionnées en fonction des intérêts des contributeurs scientifiques. Cela devrait inclure l'examen de la pertinence pour la communauté scientifique.
Les chaînes SMILES peuvent être abstraites dans une variété de méthodes :

-   Pour les molécules simples, une représentation des SMILES peut être directement traduite à l'aide de visuels.
    inspection. Ceci est généralement approprié pour les composés au début d'une liste rapportée qui contiennent les anneaux de dénominateur les plus courants.

-   Pour les molécules complexes, l'image peut être redessinée dans la version gratuite de ChemDraw, puis traduite en SMILES.

-   Pour les sources où les SMILES sont écrits et l'IUPAC n'est pas connu, les SMILES sont traduits en ChemDraw et le nom récupéré.
    Notez que certains des noms peuvent être modifiés en fonction de l'inspection humaine en faveur des noms préférés.

-   Pour les papiers polymères, les points de site ont été omis du nom et une partie de la nomenclature a été ajustée pour les noms préférés
    sur traditionnel. Par exemple : "yl" pour marquer les points de site pour les connexions polymères a été supprimé au profit d'une complexité anglaise réduite.

-   Dans le cas des radicaux, certains SMILES ont été ajustés pour supprimer la caractéristique chimique radicale car ils servent de points de connexion. Cependant, dans certains cas, la composante radicale a été maintenue, en particulier dans le cas des substituants communs du livre bleu IUPAC.

-   Les chaînes SMARTS ont été adaptées de SMILES à l'aide de RDKit (4).

* * *

# Citation

C'est en route

# Licence

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)
