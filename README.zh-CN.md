<h1 align="center">Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities </h1>

Global Chem 是一个开源图形数据库和 api，用于使用 IUPAC 作为输入和 SMILES/SMARTS 作为输出的常见和稀有化学品列表。作为
当我搜索化学无穷大时，我最需要的是我自己。

我发现这些写在历史上的列表很有用，它们来自各种不同的领域，但都是汇总的
进入有机化学家最常见的格式 (IUPAC) 和化学信息学家的通用语言 (SMILES)
模式匹配（SMARTS）。

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

# 概述

## 环球化学

#### - 入门

| [![Documentation](https://img.shields.io/badge/GitBook-Docu-lightblue)](https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/) | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1HaDvAYpaX0_2Eqx1NSl5uJTRkE_0na-L?usp=sharing) | [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.chemicalgraphtheory.com) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) |
| -------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |

#### - 验证

| [![saythanks](https://img.shields.io/badge/RDKit-100%25-fg49b4.svg)](https://www.rdkit.org/) | [![saythanks](https://img.shields.io/badge/PartialSMILES-85.7%25-fg49b4.svg)](https://github.com/baoilleach/partialsmiles) | [![saythanks](https://img.shields.io/badge/DeepSMILES-99.25%25-lm89b5.svg)](https://github.com/baoilleach/deepsmiles) | [![saythanks](https://img.shields.io/badge/SELFIES-100%25-lm89b5.svg)](https://github.com/aspuru-guzik-group/selfies) | [![saythanks](https://img.shields.io/badge/MolVS-98.5%25-lm89b5.svg)](https://github.com/mcs07/MolVS) | [![saythanks](https://img.shields.io/badge/PySMILES-99.8%25-fg49b4.svg)](https://github.com/pckroon/pysmiles) |
| -------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |

#### - 统计数据

| [![Downloads](https://pepy.tech/badge/global-chem)](https://pepy.tech/project/global-chem) | ![visitor badge](https://visitor-badge.glitch.me/badge?page_id=sulstice.global-chem) | [![Man Hours](https://img.shields.io/endpoint?url=https%3A%2F%2Fmh.jessemillar.com%2Fhours%3Frepo%3Dhttps%3A%2F%2Fgithub.com%2FSulstice%2Fglobal-chem)](https://jessemillar.com/r/man-hours) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) |
| ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |

#### - Github 操作

| [![Test System](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml) | [![pre-commit.ci status](https://results.pre-commit.ci/badge/github/Sulstice/global-chem/master.svg)](https://results.pre-commit.ci/latest/github/Sulstice/global-chem/master) | [![publish](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml) |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |

#### - 构建信息

| [![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem) | [![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) | [![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250) | [![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_shield) | [![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/) | [![Maturity badge - level 2](https://img.shields.io/badge/Maturity-Level%202%20--%20First%20Release-yellowgreen.svg)](https://github.com/tophat/getting-started/blob/master/scorecard.md) |
| ------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

#### - GlobalChemExtensions

| [![PyPI version](https://badge.fury.io/py/global-chem-extensions.svg)](https://badge.fury.io/py/global-chem-extensions) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) | [![Downloads](https://pepy.tech/badge/global-chem-extensions)](https://pepy.tech/project/global-chem-extensions) |
| ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- |

### 扩展演示卷轴

| 应用                     | 介绍                                                                                                                                                         | 高级用法                    |                                                                                                                                                            |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 力场                     | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1hW0K6V5zPDFdZvFkYarbOr9wRoof2n4s?usp=sharing) | CGenFF 分子加载器和原子类型相似性    | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1OEm1dACd_wkw_JQITemy18pgfKdE4XlX?usp=sharing) |
| 生物信息学                  | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sCNw9FIcFfTEgmrdokADXgJGnsPEOGbX?usp=sharing) | 使用 Bostrom 算法按 PDB 过滤配体 | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1a3Ys1rpqFzxBz95DQwJVHfunBxpywP7f?usp=sharing) |
| 化学信息学                  | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1z0ilrakoRJ8maapMNHwtPf83pKK1THST?usp=sharing) |                         |                                                                                                                                                            |
| 量子化学                   | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1BGLQphP1IMLndeyavHM_6_qXQJI_7gFU?usp=sharing) | 绘制量子理论和基集与小分子的哈密顿量      | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/11JGsen912TmyR5Dds-Opon5cRkrtVoT7?usp=sharing) |
| development_operations | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1v_yXPXPilbWUZGkel_yekBr3EFnIUNUL?usp=sharing) |                         |                                                                                                                                                            |

# 安装

GlobalChem 将通过 PyPi 分发，随着树及其扩展的增长，我们可以将其扩展到其他软件
无论您使用什么，都可以访问它。或者，您可以浏览一下源代码并复制/粘贴
它自己。

```bash

pip install global-chem

```

如果您想安装扩展包以获得额外的功能，每个应用程序可以相互独立安装，或者您可以使用`all`:

```bash

pip install 'global-chem[graphing]'
pip install 'global-chem[forcefields]'
pip install 'global-chem[bioinformatics]'
pip install 'global-chem[cheminformaitcs]'
pip install 'global-chem[quantum_chemistry]'
pip install 'global-chem[development_operations]'
pip install 'global-chem[all]'

```

# 快速开始

这里我们加载`global-chem[cheminformatics]`扩展包和`GlobalChem`树。我们从畅销书 Pihkal 中提取 SMILES，并对化学列表进行化学信息学主成分分析。

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

# 环球化学

### 规则

Graph Network (GN) 附带了一些规则，这些规则目前使开发人员的软件工程更容易。

-   必须有根节点。
-   添加节点时，必须连接每个节点。
-   要删除一个节点，它不能有任何子节点。

深度图网络 (DGN) 还附带了一些规则以使实施更容易：

-   必须有一个 1 的根节点标记为您的“输入”节点。
-   添加层时，所有节点都将作为子层添加到所有先前的层中。 （Folk 可以使用移除节点功能来执行 dropouts）。

### 知识图谱

只要没有依赖关系，初始化类就可以了！世界上所有常见和稀有的群体
由你处置。

<p align="center">
  <img width="850" height="650" src="images/figures/figure_1_new.png">
</p>

### 节点贡献者

如果您想选择自己或其他人，请遵循节点贡献指南。

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

| 化学品清单                 | 条目数 | 参考                                                                                                                                                                                                                            |
| --------------------- | --- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 氨基酸                   | 20  | 常识                                                                                                                                                                                                                            |
| 必需维生素                 | 13  | 常识                                                                                                                                                                                                                            |
| 常用有机溶剂                | 42  | 富尔默，格雷戈里 R.，等人。 “痕量杂质的 NMR 化学位移：与有机金属化学家相关的氘化溶剂中的常见实验室溶剂、有机物和气体。”有机金属，卷。 29，没有。 9，2010 年 5 月，第 2176-79 页。                                                                                                                     |
| 打开微笑                  | 94  | OpenSMILES 主页。[HTTP://opens Miles.org/](http://opensmiles.org/).                                                                                                                                                              |
| IUPAC 蓝皮书（CRC 手册）2003 | 333 | 化学橡胶公司。 CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edited by David R. Lide, 85. ed, CRC Press, 2004。                                                                         |
| 毒品戒指                  | 92  | 泰勒，理查德 D.，等人。 “毒品戒指。”药物化学杂志，第一卷。 57，没有。 14，2014 年 7 月，第 5845-59 页。 ACS出版物，[HTTPS://do i.org/10.1021/聚美4017625](https://doi.org/10.1021/jm4017625).                                                                            |
| 阶段 2 杂环               | 19  | 布劳顿、霍华德 B. 和伊恩 A. 沃森。 “用于药物设计的杂环的选择。”分子图形与建模杂志，第一卷。 23，没有。 1，2004 年 9 月，第 51-58 页。考研，[HTTPS://do i.org/10.1016/就.贱买贵卖.2004.03.016](https://doi.org/10.1016/j.jmgm.2004.03.016).                                               |
| 特权脚手架                 | 47  | 韦尔施，马修 E.，等人。 “图书馆设计和药物发现的特权支架。”化学生物学的当前观点，第一卷。 14，没有。 3，2010 年 6 月，第 347-61 页。PubMed，[HTTPS://do i.org/10.1016/就.才不怕.2010.02.018](https://doi.org/10.1016/j.cbpa.2010.02.018).                                               |
| 普通弹头                  | 29  | Gehringer、Matthias 和 Stefan A. Laufer。 “用于靶向共价抑制剂的新兴和重新出现的弹头：在药物化学和化学生物学中的应用。”药物化学杂志，卷。 62，没有。 12，2019 年 6 月，第 5673-724 页。 ACS出版物，[HTTPS://do i.org/10.1021/ACS.就么的车门.8不01153](https://doi.org/10.1021/acs.jmedchem.8b01153). |
| 常见的聚合物重复单元            | 78  | Hiorns, R. C. 等人。 “聚合物命名法简要指南（IUPAC 技术报告）。”纯化学和应用化学，卷。 84，没有。 10，2012 年 10 月，第 2167-69 页。[HTTPS://do i.org/10.1351/PAC-rep-12-03-05](https://doi.org/10.1351/PAC-REP-12-03-05).                                               |
| 常见的 R 组替换             | 499 | 竹内，康介等人。 “药物化学的 R 组替代数据库。”未来科学 OA 卷。 7，没有。 2021 年 9 月 8 日，第 8 页。 FSO742。未来科学网（Atypon），[HTTPS://do i.org/10.2144/风骚哦啊-2021-0062](https://doi.org/10.2144/fsoa-2021-0062).                                                      |
| 激酶的亲电弹头               | 24  | Petri、László 等人。 “用于绘制蛋白质激酶中易处理半胱氨酸的反应性和可及性的亲电弹头库。”欧洲药物化学杂志，第一卷。 207，2020 年 12 月，第112836.考研，[HTTPS://do i.org/10.1016/就.阿胶么吃.2020.112836](https://doi.org/10.1016/j.ejmech.2020.112836).                                      |
| 激酶的特权支架               | 29  | 胡华斌等。 “竞争性和变构激酶抑制剂的系统比较揭示了共同的结构特征。”欧洲药物化学杂志，第一卷。 214，2021 年 3 月，第113206. 科学直接,[HTTPS://do i.org/10.1016/就.阿胶么吃.2021.113206](https://doi.org/10.1016/j.ejmech.2021.113206).                                                    |
| BRaf 抑制剂              | 54  | Agianian、Bogos 和 Evripidis Gavathiotis。 “BRAF 抑制剂在癌症中的最新见解。”药物化学杂志，第一卷。 61，没有。 14，2018 年 7 月，第 5775-93 页。 ACS出版物，[HTTPS://do i.org/10.1021/ACS.就么的车门.7不01306](https://doi.org/10.1021/acs.jmedchem.7b01306).                  |
| 常见氨基酸保护基团             | 346 | Isidro-Llobet，阿尔伯特等人。 “氨基酸保护基团。”化学评论，第一卷。 109，没有。 6，2009 年 6 月，第 2455-504 页。 DOI.org（交叉引用），[HTTPS://do i.org/10.1021/成人800323是](https://doi.org/10.1021/cr800323s).                                                           |
| 新兴的全氟烷基               | 27  | 佩尔奇，凯瑟琳 E.，等人。 “PFAS 健康影响数据库：系统证据图协议。”环境国际，第一卷。 130，2019 年 9 月，第104851. 科学直接,[HTTPS://do i.org/10.1016/就.恶女int.2019.05.045](https://doi.org/10.1016/j.envint.2019.05.045).                                                    |
| 粘土吸附用化学品              | 33  | 奥尔，明日香 A.，等人。 “结合实验等温线、简约模拟和模型来理解和预测蒙脱石粘土上的化学吸附。” ACS 欧米茄，第一卷。 6，没有。 22，2021 年 6 月，第 14090-103 页。考研，[HTTPS://do i.org/10.1021/ACS omega.1从00481](https://doi.org/10.1021/acsomega.1c00481).                                   |
| 大麻素                   | 63  | 特纳，卡尔顿 E.，等人。 “大麻苜蓿 L. XVII 的成分。天然成分回顾。”天然产物杂志，第一卷。 43，没有。 2，1980 年 3 月，第 169-234 页。 ACS出版物，[HTTPS://do i.org/10.1021/哪怕50008啊001](https://doi.org/10.1021/np50008a001).                                                      |
| 附表 1 美国麻醉品            | 240 | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                                                                 |
| 附表 2 美国麻醉品            | 60  | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                                                                 |
| 附表 3 美国麻醉品            | 22  | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                                                                 |
| 附表 4 美国麻醉品            | 77  | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                                                                 |
| 附表 5 美国麻醉品            | 8   | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                                                                 |
| 皮卡尔                   | 179 | 舒尔金、亚历山大 T. 和安舒尔金。 Pihkal：化学爱情故事。 1. ed.，8. 印刷，Transform，2010。                                                                                                                                                                |
| 辅料西咪替丁和阿昔洛韦           | 14  | Vaithianathan、Soundarya 等人。 “常见赋形剂对生物药剂学分类系统第 3 类药物西咪替丁和阿昔洛韦口服药物吸收的影响。”药学杂志，第一卷。 105，没有。 2，2016 年 2 月，第 996-1005 页。考研，[HTTPS://do i.org/10.1002/吉普赛.24643](https://doi.org/10.1002/jps.24643).                                |
| 镍双膦酸盐配体               | 不适用 | Clevenger、Andrew L. 等人。 “在镍催化中使用二齿膦作为配体的趋势。”化学评论，第一卷。 120，没有。 13，2020 年 7 月，第 6124-96 页。 DOI.org（交叉引用），[HTTPS://do i.org/10.1021/ACS.Chem Rev.9不00682](https://doi.org/10.1021/acs.chemrev.9b00682).                          |
| 常见的正则表达式模式            | 1   |                                                                                                                                                                                                                               |

</details>                            

# GlobalChemExtensions

您可以使用多种工具来浏览和分析数据，并且可以在 google colab 演示或 Gitbook 文档中找到不同应用程序的完整列表。使用 plotly 和 bokeh 设计的数据可视化扩展的演示如下所示：

<p align="center">
  <img width="800" height="600" src="https://raw.githubusercontent.com/Sulstice/global-chem/master/images/figures/figure_10.png">
</p>

<details><summary><h3>Extension List</h1><br/></summary>

| 延期              | 描述                                                                            | 应用                     |
| --------------- | ----------------------------------------------------------------------------- | ---------------------- |
| GlobalChem 化学实体 | GlobalChem 具有内部分子对象，具有关联的所有常见属性并转换为 SMILES                                    | 力场                     |
| GlobalChem 生物实体 | GlobalChem 具有内部 DNA/RNA/蛋白质/分子对象，具有关联的所有常见属性并转换为 SMILES                       | 生物信息学                  |
| 可视化 DNA/RNA 链   | 可视化 DNA 和 RNA 链并为其添加标签                                                        | 生物信息学                  |
| 力场分子            | GlobalChem 可以将 CGenFF 和 GaFF2 文件作为对象解析、操作和写入                                  | 力场                     |
| PDF 生成和解析       | GlobalChem 可以将 SMILES 生成 PDF 并将 PDF 转换为 SMILES                                | 化学信息学                  |
| 微笑验证            | GlobalChem 与 PySMILES、DeepSMILES、PartialSmiles、SELFIES、MolVS 连接，用于验证 SMILES 集 | 化学信息学                  |
| 微笑质子状态          | GlobalChem 可以采用一组化合物并预测 SMILES 串在一定 pH 值范围内的质子化状态                             | 化学信息学                  |
| 开源数据库监控         | GlobalChem 使用 Uptime-Cheminformatics 跟踪开源化学数据                                 | development_operations |
| Networkx 软件适配器  | GlobalChem Network 可以转换为 NetworkX Graph Objects                               | 化学信息学                  |
| SMARTS 模式验证     | GlobalChem 使用 MiniFrag 数据库测试 SMARTS 字符串在功能组选择中的准确性                            | 化学信息学                  |
| 主成分分析           | GlobalChem 可以轻松解释 SMILES、指纹、聚类和应用 PCA 分析，用户可以调整参数                             | 化学信息学                  |
| 药物设计过滤器         | GlobalChem 可以根据通用药物设计过滤规则过滤化合物                                                | 化学信息学                  |
| 深层散射分析          | 为了可视化分子组之间的关系，GlobalChem 提供了平行坐标图生成                                           | 化学信息学                  |
| 旭日形径向分析         | GlobalChem 提供了一种旭日机制，允许用户观察化合物组与公共组的关系                                        | 化学信息学                  |
| 图形模板            | GlobalChem 提供图形模板以帮助更快地进行数据分析，目前唯一提供的是 Plotly                                 | 化学信息学                  |
| CGenFF 差异分数     | GlobalChem 可以根据原子类型提供两种分子之间的差异                                                | 力场                     |
| 一种热编码           | GlobalChem 拥有自己的一个基于机器学习通用列表的热门编码器和解码器                                        | 化学信息学                  |
| SMARTS 模式标识符    | GlobalChem 连接到 SMARTS Plus，可以提供不同 SMARTS 组件的可视化                               | 化学信息学                  |
| Psi4 解析器        | 提供 Psi4 输出文件的解析和提取值                                                           | 量子化学                   |
| 坐标商店            | 用于在 xyz 和 z 矩阵中分布的小分子坐标的仓库                                                    | 量子化学                   |
| 可视化分子轨道         | 可视化来自 Psi4 输出 cubeprop 的多维数据集文件                                               | 量子化学                   |

</details>

# 开源软件合规性

`GlobalChem`遵循联邦法规第 21 篇第 11 部分概述的相同原则；电子记录，
电子签名（21 CFR 第 11 部分）指导文件。由于没有关于如何处理开源软件的正式指南，我们
尝试完成要求。 FDA 认为第 11 部分适用于以下电子记录标准以及如何`GlobalChem`完成每个组件：

-   **合理性：**`GlobalChem`建立在使用阅读和重绘从书籍和论文中提取的数据之上。它添加了一个组件
    IUPAC/SMILES/SMARTS 字符串以电子方式存储它，从而为它的数据提供独特的组件。记录是开源的
    并由存储库的维护者和开源社区反馈进行适当的版本控制。`GlobalChem`的目的仍然未知，因为它进入了开源部署。我们已经构建了扩展功能
    一个单独的包`GlobalChemExtensions`这确实取决于`GlobalChem`.由于每个版本都被适当地打包，
    需要依赖一个版本，然后它的软件才可用`Github`和`PyPi`.标准操作程序 (SOP)
    可以从维护的扩展实用程序文档中提交`Gitbook`.

-   **验证：**`GlobalChem`遵循良好自动化制造规范 (GAMP) 类别 3，即“安装后使用的软件”
    并且可能是“可配置的”。`GlobalChem`测试来自内部，文档作为最终测试
    功能，因为这是用户测试最多的东西，因为我们依赖开源。持续集成 (CI)
    系统还同时构建以作为基本功能测试`GlobalChem`图网络。存储的数据
    由该领域的专家维护，但如果发现错误，可能会根据社区反馈进行更改。

-   **审计追踪：**`GlobalChem`受版本控制`Git`并托管在微软的平台上`Github`.`GlobalChem`遵循语义
    模式的版本控制`X1.X2.X3`:`X1`用测试、文档和平均值标记正式的稳定版本
    对软件或功能进行大的重构，`X2`意味着添加了一个新功能，无论是否有测试和文档，但
    如此迭代。`X3`表示要添加到函数中的“热”修复（一个简单的错误）、小功能或附加参数
    ，或迭代数据。

-   **旧系统：**`GlobalChem`自第一个版本发布以来已经运行了近 2 年`0.3.0`2020 年 5 月。`GlobalChem`在开源社区中建立了完整的路径，每个版本都被编目并且对所有人可见。这满足
    规则概述了确定遗留系统。我们使用社交媒体平台（Twitter、Github、LinkedIn）提供的社区反馈
    作为书面证据和理由`GlobalChem`适合化学信息学的预期用途。

-   **记录副本：**`GlobalChem`有记录存储在`Github`对于可以导出为各种格式的软件
    微软。对于文档，它托管在`Gitbook`并根据软件进行版本控制。每一个“书”
    可以导出为适用于 FDA 提交的便携式数据格式 (PDF)。

-   **记录保留：**`GlobalChem`具有控制为唯一 ID (UUID) 的文档版本记录，该唯一 ID (UUID) 用作其标识符
    对于存储在`Gitbook`.每个版本都存储为降价文件，并在需要时转换为 PDF。

`GlobalChem`具有 Mozilla 公共许可证 2.0 版。`GlobalChem`允许您在更大的工作中使用该软件，并且
如果您愿意，可以通过修改对其进行扩展。意外情况是，如果您安装`GlobalChem`并发布新软件
那么您必须遵循我们的开源社区许可证中安装的相同原则。

# 数据采集

参考文献和相关化合物列表是根据科学贡献者的兴趣选择的。这应包括考虑与科学界的相关性。
SMILES 字符串可以用多种方法进行抽象：

-   对于简单的分子，可以使用视觉直接翻译 SMILES 的一种表示
    检查。这通常适用于报告列表开头包含最常见分母环的化合物。

-   对于复杂的分子，图像可以在免费版本的 ChemDraw 中重新绘制，然后翻译成 SMILES。

-   对于写入 SMILES 且不知道 IUPAC 的来源，SMILES 将被翻译成 ChemDraw 并检索名称。
    请注意，某些名称可能会根据人工检查进行修改，以支持首选名称。

-   对于聚合物纸，名称中省略了站点点，并且根据首选名称调整了一些命名法
    超过传统。例如：用于标记聚合物连接点的“yl”已被删除，以降低英语的复杂性。

-   在自由基的情况下，一些 SMILES 被调整以去除自由基化学特征，因为它们用作连接点。然而，在某些情况下，自由基成分得以保留，特别是在 IUPAC 蓝皮书常见取代基的情况下。

-   SMARTS 字符串是使用 RDKit (4) 从 SMILES 改编而来的。

* * *

# 引文

它正在路上

# 许可

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)
