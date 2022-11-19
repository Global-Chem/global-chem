<h1 align="center">Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities </h1>

Global Chemは入力にはIUPACを使い、SMILES/SMARTSを出力する、一般的な、または希少化学物資リストのためのオープンソースグラフデータベース、及びそのためのAPIです。膨大な化学物質を検索する際に利用できます。

私は、これらの異なる分野から集められた、有機化学者の共通フォーマットであるIUPACや、ケモインフォマティシャンの共通言語であるSMILESに集約されているこれらのリストがとても有用であるということに気がついたのです。

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

オーバービュー
=============

## GlobalChem

#### - Getting Started

|[![Documentation](https://img.shields.io/badge/GitBook-Docu-lightblue)](https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/) | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1HaDvAYpaX0_2Eqx1NSl5uJTRkE_0na-L?usp=sharing) | [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.chemicalgraphtheory.com)| [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) |
|-|-|-|-|

#### - Validation 

|[![saythanks](https://img.shields.io/badge/RDKit-100%25-fg49b4.svg)](https://www.rdkit.org/)| [![saythanks](https://img.shields.io/badge/PartialSMILES-85.7%25-fg49b4.svg)](https://github.com/baoilleach/partialsmiles) | [![saythanks](https://img.shields.io/badge/DeepSMILES-99.25%25-lm89b5.svg)](https://github.com/baoilleach/deepsmiles) | [![saythanks](https://img.shields.io/badge/SELFIES-100%25-lm89b5.svg)](https://github.com/aspuru-guzik-group/selfies) | [![saythanks](https://img.shields.io/badge/MolVS-98.5%25-lm89b5.svg)](https://github.com/mcs07/MolVS) | [![saythanks](https://img.shields.io/badge/PySMILES-99.8%25-fg49b4.svg)](https://github.com/pckroon/pysmiles) |
|-|-|-|-|-|-|

#### - Statistics

| [![Downloads](https://pepy.tech/badge/global-chem)](https://pepy.tech/project/global-chem) |![visitor badge](https://visitor-badge.glitch.me/badge?page_id=sulstice.global-chem) | [![Man Hours](https://img.shields.io/endpoint?url=https%3A%2F%2Fmh.jessemillar.com%2Fhours%3Frepo%3Dhttps%3A%2F%2Fgithub.com%2FSulstice%2Fglobal-chem)](https://jessemillar.com/r/man-hours) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)|
|-|-|-|-|

#### - Github Actions

|[![Test System](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml)|[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/Sulstice/global-chem/master.svg)](https://results.pre-commit.ci/latest/github/Sulstice/global-chem/master)| [![publish](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml) | [![Translate README](https://github.com/Sulstice/global-chem/actions/workflows/translate_readme.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/translate_readme.yml) |
|-|-|-|-|

#### - Build Information

| [![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem) | [![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)| [![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250) | [![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_shield) | [![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/) | [![Maturity badge - level 2](https://img.shields.io/badge/Maturity-Level%202%20--%20First%20Release-yellowgreen.svg)](https://github.com/tophat/getting-started/blob/master/scorecard.md) |
|-|-|-|-|-|-|-|

#### - GlobalChemExtensions

| [![PyPI version](https://badge.fury.io/py/global-chem-extensions.svg)](https://badge.fury.io/py/global-chem-extensions) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) | [![Downloads](https://pepy.tech/badge/global-chem-extensions)](https://pepy.tech/project/global-chem-extensions) |
|-|-|-|

### Extension Demo Reel

| Application              | Introduction | Advanced Usage  | |
| ---------------------- | ---- | ---------------------- |-|
| forcefields            | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1hW0K6V5zPDFdZvFkYarbOr9wRoof2n4s?usp=sharing)  | CGenFF Molecule Loader and Atom Type Similarity | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1OEm1dACd_wkw_JQITemy18pgfKdE4XlX?usp=sharing) |
| bioinformatics         | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sCNw9FIcFfTEgmrdokADXgJGnsPEOGbX?usp=sharing)  | Use the Bostrom Algorithm to Filter Ligands By PDB  | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1a3Ys1rpqFzxBz95DQwJVHfunBxpywP7f?usp=sharing) |
| cheminformatics        | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1z0ilrakoRJ8maapMNHwtPf83pKK1THST?usp=sharing)   |            | |
| quantum_chemistry      | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1BGLQphP1IMLndeyavHM_6_qXQJI_7gFU?usp=sharing) | Plot Quantum Theory and Basis Set versuses the Hamiltonian of small molecules | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/11JGsen912TmyR5Dds-Opon5cRkrtVoT7?usp=sharing) |
| development_operations | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1v_yXPXPilbWUZGkel_yekBr3EFnIUNUL?usp=sharing) |

インストール 
==========

GlobalChemはPyPiから配布されます。そしてパッケージとその追加機能が更新されるに連れて様々なソフトからアクセスできるように拡張できるようになるでしょう。またはソースコードを見てコピー＆ペーストして使うこともできます。

```bash

pip install global-chem

```

もし追加機能をインストールしたい場合は、それぞれのアプリケーションを個別にインストールすることもできますし、`all`:　をつけて全部を一気にインストールすることもできます。

```bash

pip install 'global-chem[graphing]'
pip install 'global-chem[forcefields]'
pip install 'global-chem[bioinformatics]'
pip install 'global-chem[cheminformaitcs]'
pip install 'global-chem[quantum_chemistry]'
pip install 'global-chem[development_operations]'
pip install 'global-chem[all]'

```

クイックスタート
==============

さあ、まずは`global-chem[cheminformatics]`拡張機能と、`GlobalChem`をロードします。そして有名な書籍PihkalからSMILESを抜き出し、これらのリストについてケモインフォマティクス的な主成分分析を実施します。

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

GlobalChem
==========

### ルール


グラフネットワーク（GN)は開発者のソフトウェアエンジニアリングを容易にするためのいくつかのルールに基づいて構成されます。

- ルートノードは必須である。
- ノードを追加する場合、すべてのノードは結合されていなければならない。 
- ノードを削除する場合、そのノードが子ノードを持っていてはならない。

Deepグラフネットワーク(DGN)も同様に、実装を容易にするためのいくつかのルールに基づいて構成されます。

- 入力ノードであると示すための１がルートノードに必要。 
- 層を追加するときは、すべてのノードが前層に子として追加される。（ノードの特徴はドロップアウト実施することで取り除くことができます。）

### Knowledge グラフ


依存関係はありません。クラスを初期化すれば、それで完了! 世界の共通グループと希少グループのKnowledge グラフを扱えるようになります。

<p align="center">
  <img width="850" height="650" src="images/figures/figure_1_new.png">
</p>

### Nodes 貢献者


自薦他薦を問わず、以下のノード貢献に関するガイドラインを確認し投稿して下さい。

```
'global_chem': Node,
'emerging_perfluoroalkyls': EmergingPerFluoroAlkyls,                      # Asuka Orr & Suliman Sharif
'montmorillonite_adsorption': MontmorilloniteAdsorption,      `GlobalChem`nhibitors,  # Shaoqi Zhao & Suliman Sharif
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
```

<details><summary><h3>Node リスト</h1><br/></summary>

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
| Cannabinoids                        | 63           | Turner, Carlton E., et al. “Constituents of Cannabis Sativa L. XVII. A Review of the Natural Constituents.” Journal of Natural Products, vol. 43, no. 2, Mar. 1980, pp. 169–234. ACS Publications, https://doi.org/10.1021/np50008a001.                                                                              |
| Schedule 1 United States Narcotics  | 240          | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 2 United States Narcotics  | 60           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 3 United States Narcotics  | 22           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 4 United States Narcotics  | 77           | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Schedule 5 United States Narcotics  | 8            | ECFR :: 21 CFR Part 1308 - Schedules.                                                                                                                                                                                                                                                                                |
| Pihkal                              | 179          | Shulgin, Alexander T., and Ann Shulgin. Pihkal: A Chemical Love Story. 1. ed., 8. print, Transform, 2010.                                                                                                                                                                                                            |
| Excipients Cimetidine & Acyclovir   | 14           | Vaithianathan, Soundarya, et al. “Effect of Common Excipients on the Oral Drug Absorption of Biopharmaceutics Classification System Class 3 Drugs Cimetidine and Acyclovir.” Journal of Pharmaceutical Sciences, vol. 105, no. 2, Feb. 2016, pp. 996–1005. PubMed, https://doi.org/10.1002/jps.24643.                |
| Nickel Bidendate Phosphine Ligands  | N/A          | Clevenger, Andrew L., et al. “Trends in the Usage of Bidentate Phosphines as Ligands in Nickel Catalysis.” Chemical Reviews, vol. 120, no. 13, July 2020, pp. 6124–96. DOI.org (Crossref), https://doi.org/10.1021/acs.chemrev.9b00682.                                                                              |
| Common Regex Patterns               | 1            |                                                                                                     

</details>                            

GlobalChemExtensions
====================

google colab demoやGitbookドキュメンテーションから、様々なアプリケーションのリストを利用し他データの可視化や、解析のための様々なツールが使えることがわかるでしょう。以下のデータ可視化のための拡張機能はplotly, bokehを利用してデザインされています。:

<p align="center">
  <img width="800" height="600" src="https://raw.githubusercontent.com/Sulstice/global-chem/master/images/figures/figure_10.png">
</p>

<details><summary><h3>Extension List</h1><br/></summary>

| Extension                       | Description                                                                                                             | Appplication   |
|---------------------------------|-------------------------------------------------------------------------------------------------------------------------|------------------|
| GlobalChem Chemical Entities    | GlobalChem has internal Molecule objects with all common attributes associated and conversion to SMILES                 | forcefields       |
| GlobalChem Biological Entities  | GlobalChem has internal DNA/RNA/Protein/Molecule objects with all common attributes associated and conversion to SMILES | bioinformatics   |
  | Visualize DNA/RNA Strands  | Visualize DNA and RNA Strands and add labels to them | bioinformatics   |
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

</details>

オープンソースソフトウェアコンプライアンス
===================================


`GlobalChem`は、連邦規則集のタイトル21のパート11電子署名（21 CFR Part 11）ガイダンス文書に概説されている原則に従います。オープンソースソフトウェアがどのように扱われるべきかについての公式なガイドラインがないため、私達は要件を満たすことこを試みています。FDAは、電子記録についてPart11の基準を適用することを検討しており、
`GlobalChem`はそれぞれのコンポーネントでどのように達成しようとしているかというと:


- **Plausabilitiy/妥当性:** `GlobalChem`は書籍や論文から読み込んだり、再描画したデータによって構成されています。IUPAC/SMILES/SMARTS文字列は電子的に保存し、データに独自の要素を持たせています。レコードはオープンソースであり、コミュニティからのFBなども含めてメインテナーによりリポジトリのバージョン管理がなされます。
`GlobalChem`のオープンソースの展開としての目的はまだわかりません。`GlobalChem`に依存する別のパッケージ`GlobalChemExtensions`を開発しています。各バージョンを適切にパッケージングしているため、依存するソフトが必要な場合は`Github`や`PyPi`から適切なバージョンのソフトが利用可能です。標準手順書（SOP)は`Gitbook`で管理しているドキュメントに投稿してあります。


- **Validation/検証:** `GlobalChem`は自動製造実践規範（GAMP）カテゴリ3、"インストールしたまま利用する設定可能なソフトウェア"に従います。`GlobalChem`のテストはドキュメントが究極のテストとして機能します。なぜなら、私達はオープンソースに依存しており、ユーザーが一番テストを実施するからです。継続的インテグレーション（CI)システムも構築し、`GlobalChem`グラフネットワークの基本機能のテストとして提供しています。保存されているデータはこの分野のエキスパートがメンテナンスしていますが、コミュニティーからのエラーに関するフィードバックによって変更する場合があります。


- **Audit Trail/監査証跡:** `GlobalChem` は`Git`と、マイクロソフト社によってホストされている`Github`を利用してバージョン管理をしています。 `GlobalChem` はセマンティックバージョニングのスキーマに従います。 `X1.X2.X3`:の `X1`は公式のテストとドキュメントを含む安定版リリースであること、大幅なソフトや機能のリファクタリングしたものであることを示します。 `X2` は新機能を追加したことを示しますがｍドキュメントやテストはある場合とない場合があります。 `X3`は"hot"fix（何かの軽微なバグ対応）を示し、小さな機能やパラメータの追加、データの追加などが行われます。


- **Legacy Systems/レガシーシステム:** `GlobalChem` はおよそ2年前に最初のバージョン`0.3.0`が2020年5月にリリースされてからずっと稼働しています。`GlobalChem`はオープンソースコミュニティによって構築されており、各バージョンはカタログ化されすべて見えるようになっています。これはレガシーシステムを決定するためのルールの概要を満たしています。私達はSNS（Twitter, Githb, LinkdeIn)を通じたコミュニティーからのフィードバックを`GlobalChem`がケモインフォマティクスの用途に適したものであるということの文書化された証拠、判断のために利用します。


- **Copies of Recordsレコードのコピー:** `GlobalChem`は`Github`上に保存されたレコードを有しておりMicrosoftから提供される様々なフォーマットとしてエクスポート可能です。文書化は、`Gitbook`にてホストされ、ソフトに合わせてバージョンコントロールしています。それぞれの"book"はFDA提出用のPDF形式でエクスポートすることができます。


- **Record Retention/レコードの保持:** `GlobalChem` は文書化されたレコードを有しており、それらは一意なID（UUID）を識別子として`Gitbook`にてバージョン管理されています。各バージョンはマークダウンファイルで保存されており必要に応じてPDF形式に変換可能です。 


`GlobalChem`はMozilla Public License version 2.0にて利用できます。`GlobalChem`は貴方のより大きなソフトの中で利用することや、機能を拡張することも許諾しています。もしあなたが`GlobalChem`をインストールして、新しいソフトウェアをリリースする場合、オープンソースコミュニティのために私たちのライセンスにインストールされているのと同じ原則に従わなければならないという不測の事態が発生する可能性があります。オープンソースコミュニティのためのライセンスにインストールされているのと同じ原則に従わなければなりません。


データコレクション
===============

リファレンス、そして関連する化合物リストはサイエンティフィックな貢献者らの興味に基づいて選ばれています。SMILES文字列は様々な方法で抽象化されているかもしれません:

- 単分な分子におけるSMILES表現はビジュアルインスペクションによって直接変換できます。これは一般的に、最も一般的な共通環構造を含むリストの最初の方にある化合物に適しています。

- 複雑な分子の画像は無料版のChemDrawで書き直し、それをSMILESに変換しました。

- SMILESが書かれているがIUPAC名が知られていないソースでは、ChemDrawを使って化合物名に変換しました。何個かの名前は、人に好まれる形に修飾しています。

- 高分子の論文については、サイト点を名前から省略し、一部の命名法を伝統的に好まれる名称に調整しています。例えばサイト点の印となる'yl'は英語の複雑性を減らすために取り除いています。

- ラジカルを扱う場合、いくつかのSMILESではラジカルを取り除いています。なぜならこれは結合点を意味してしまうからです。しかし、特にIUPACブルーブックの一般的な置換基の場合，ラジカル成分が維持されているケースもあった。

- SMARTS文字列はRDKit(4)を用いてSMILESから生成しました。 

* * * * *

Citation
========

It's on it's way


Licensing
=========

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)
