<h1 align="center">Global-Chem: A Chemical Knowledge Graph of common small molecules and their IUPAC/SMILES/SMARTS for selection of compounds relevant to diverse chemical communities </h1>

Global Chem 是一個開源圖形數據庫和 api，用於使用 IUPAC 作為輸入和 SMILES/SMARTS 作為輸出的常見和稀有化學品列表。作為
當我搜索化學無窮大時，我最需要的是我自己。

我發現這些寫在歷史上的列表很有用，它們來自各種不同的領域，但都是匯總的
進入有機化學家最常見的格式 (IUPAC) 和化學信息學家的通用語言 (SMILES)
模式匹配（SMARTS）。

<p align="center">
  <img width="800" height="400" src="images/globalchemlogo.png">
</p>

# 概述

## 環球化學

#### - 入門

| [![Documentation](https://img.shields.io/badge/GitBook-Docu-lightblue)](https://sulstice.gitbook.io/globalchem-your-chemical-graph-network/) | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1HaDvAYpaX0_2Eqx1NSl5uJTRkE_0na-L?usp=sharing) | [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.chemicalgraphtheory.com) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) |
| -------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |

#### - 驗證

| [![saythanks](https://img.shields.io/badge/RDKit-100%25-fg49b4.svg)](https://www.rdkit.org/) | [![saythanks](https://img.shields.io/badge/PartialSMILES-85.7%25-fg49b4.svg)](https://github.com/baoilleach/partialsmiles) | [![saythanks](https://img.shields.io/badge/DeepSMILES-99.25%25-lm89b5.svg)](https://github.com/baoilleach/deepsmiles) | [![saythanks](https://img.shields.io/badge/SELFIES-100%25-lm89b5.svg)](https://github.com/aspuru-guzik-group/selfies) | [![saythanks](https://img.shields.io/badge/MolVS-98.5%25-lm89b5.svg)](https://github.com/mcs07/MolVS) | [![saythanks](https://img.shields.io/badge/PySMILES-99.8%25-fg49b4.svg)](https://github.com/pckroon/pysmiles) |
| -------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |

#### - 統計數據

| [![Downloads](https://pepy.tech/badge/global-chem)](https://pepy.tech/project/global-chem) | ![visitor badge](https://visitor-badge.glitch.me/badge?page_id=sulstice.global-chem) | [![Man Hours](https://img.shields.io/endpoint?url=https%3A%2F%2Fmh.jessemillar.com%2Fhours%3Frepo%3Dhttps%3A%2F%2Fgithub.com%2FSulstice%2Fglobal-chem)](https://jessemillar.com/r/man-hours) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) |
| ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |

#### - Github 操作

| [![Test System](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/continous_integration.yml) | [![pre-commit.ci status](https://results.pre-commit.ci/badge/github/Sulstice/global-chem/master.svg)](https://results.pre-commit.ci/latest/github/Sulstice/global-chem/master) | [![publish](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml/badge.svg)](https://github.com/Sulstice/global-chem/actions/workflows/publish_package.yml) |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |

#### - 構建信息

| [![PyPI version](https://badge.fury.io/py/global-chem.svg)](https://badge.fury.io/py/global-chem) | [![Coverage Status](https://coveralls.io/repos/github/Sulstice/global-chem/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/global-chem?branch=master) | ![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem) | [![DOI](https://zenodo.org/badge/259046250.svg)](https://zenodo.org/badge/latestdoi/259046250) | [![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_shield) | [![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/) | [![Maturity badge - level 2](https://img.shields.io/badge/Maturity-Level%202%20--%20First%20Release-yellowgreen.svg)](https://github.com/tophat/getting-started/blob/master/scorecard.md) |
| ------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

#### - GlobalChemExtensions

| [![PyPI version](https://badge.fury.io/py/global-chem-extensions.svg)](https://badge.fury.io/py/global-chem-extensions) | [![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0) | [![Downloads](https://pepy.tech/badge/global-chem-extensions)](https://pepy.tech/project/global-chem-extensions) |
| ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- |

### 擴展演示捲軸

| 應用                     | 介紹                                                                                                                                                         | 高級用法                    |                                                                                                                                                            |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 力場                     | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1hW0K6V5zPDFdZvFkYarbOr9wRoof2n4s?usp=sharing) | CGenFF 分子加載器和原子類型相似性    | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1OEm1dACd_wkw_JQITemy18pgfKdE4XlX?usp=sharing) |
| 生物信息學                  | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sCNw9FIcFfTEgmrdokADXgJGnsPEOGbX?usp=sharing) | 使用 Bostrom 算法按 PDB 過濾配體 | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1a3Ys1rpqFzxBz95DQwJVHfunBxpywP7f?usp=sharing) |
| 化學信息學                  | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1z0ilrakoRJ8maapMNHwtPf83pKK1THST?usp=sharing) |                         |                                                                                                                                                            |
| 量子化學                   | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1BGLQphP1IMLndeyavHM_6_qXQJI_7gFU?usp=sharing) | 繪製量子理論和基集與小分子的哈密頓量      | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/11JGsen912TmyR5Dds-Opon5cRkrtVoT7?usp=sharing) |
| development_operations | [![Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1v_yXPXPilbWUZGkel_yekBr3EFnIUNUL?usp=sharing) |                         |                                                                                                                                                            |

# 安裝

GlobalChem 將通過 PyPi 分發，隨著樹及其擴展的增長，我們可以將其擴展到其他軟件
無論您使用什麼，都可以訪問它。或者，您可以瀏覽一下源代碼並複制/粘貼
它自己。

```bash

pip install global-chem

```

如果您想安裝擴展包以獲得額外的功能，每個應用程序可以相互獨立安裝，或者您可以使用`all`:

```bash

pip install 'global-chem[graphing]'
pip install 'global-chem[forcefields]'
pip install 'global-chem[bioinformatics]'
pip install 'global-chem[cheminformaitcs]'
pip install 'global-chem[quantum_chemistry]'
pip install 'global-chem[development_operations]'
pip install 'global-chem[all]'

```

# 快速開始

這裡我們加載`global-chem[cheminformatics]`擴展包和`GlobalChem`樹。我們從暢銷書 Pihkal 中提取 SMILES，並對化學列表進行化學信息學主成分分析。

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

# 環球化學

### 規則

Graph Network (GN) 附帶了一些規則，這些規則目前使開發人員的軟件工程更容易。

-   必須有根節點。
-   添加節點時，必須連接每個節點。
-   要刪除一個節點，它不能有任何子節點。

深度圖網絡 (DGN) 還附帶了一些規則以使實施更容易：

-   必須有一個 1 的根節點標記為您的“輸入”節點。
-   添加層時，所有節點都將作為子層添加到所有先前的層中。 （Folk 可以使用移除節點功能來執行 dropouts）。

### 知識圖譜

只要沒有依賴關係，初始化類就可以了！世界上所有常見和稀有的群體
由你處置。

<p align="center">
  <img width="850" height="650" src="images/figures/figure_1_new.png">
</p>

### 節點貢獻者

如果您想選擇自己或其他人，請遵循節點貢獻指南。

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

| 化學品清單                 | 條目數 | 參考                                                                                                                                                                                |
| --------------------- | --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 氨基酸                   | 20  | 常識                                                                                                                                                                                |
| 必需維生素                 | 13  | 常識                                                                                                                                                                                |
| 常用有機溶劑                | 42  | 富爾默，格雷戈里 R.，等人。 “痕量雜質的 NMR 化學位移：與有機金屬化學家相關的氘化溶劑中的常見實驗室溶劑、有機物和氣體。”有機金屬，卷。 29，沒有。 9，2010 年 5 月，第 2176-79 頁。                                                                         |
| 打開微笑                  | 94  | OpenSMILES 主頁。<http://opensmiles.org/>.                                                                                                                                           |
| IUPAC 藍皮書（CRC 手冊）2003 | 333 | 化學橡膠公司。 CRC Handbook of Chemistry and Physics: A Ready-Reference Book of Chemical and Physical Data Edited by David R. Lide, 85. ed, CRC Press, 2004。                             |
| 毒品戒指                  | 92  | 泰勒，理查德 D.，等人。 “毒品戒指。”藥物化學雜誌，第一卷。 57，沒有。 14，2014 年 7 月，第 5845-59 頁。 ACS出版物，<https://doi.org/10.1021/jm4017625>.                                                                    |
| 階段 2 雜環               | 19  | 布勞頓、霍華德 B. 和伊恩 A. 沃森。 “用於藥物設計的雜環的選擇。”分子圖形與建模雜誌，第一卷。 23，沒有。 1，2004 年 9 月，第 51-58 頁。考研，<https://doi.org/10.1016/j.jmgm.2004.03.016>.                                                |
| 特權腳手架                 | 47  | 韋爾施，馬修 E.，等人。 “圖書館設計和藥物發現的特權支架。”化學生物學的當前觀點，第一卷。 14，沒有。 3，2010 年 6 月，第 347-61 頁。PubMed，<https://doi.org/10.1016/j.cbpa.2010.02.018>.                                               |
| 普通彈頭                  | 29  | Gehringer、Matthias 和 Stefan A. Laufer。 “用於靶向共價抑製劑的新興和重新出現的彈頭：在藥物化學和化學生物學中的應用。”藥物化學雜誌，卷。 62，沒有。 12，2019 年 6 月，第 5673-724 頁。 ACS出版物，<https://doi.org/10.1021/acs.jmedchem.8b01153>. |
| 常見的聚合物重複單元            | 78  | Hiorns, R. C. 等人。 “聚合物命名法簡要指南（IUPAC 技術報告）。”純化學和應用化學，卷。 84，沒有。 10，2012 年 10 月，第 2167-69 頁。<https://doi.org/10.1351/PAC-REP-12-03-05>.                                              |
| 常見的 R 組替換             | 499 | 竹內，康介等人。 “藥物化學的 R 組替代數據庫。”未來科學 OA 卷。 7，沒有。 2021 年 9 月 8 日，第 8 頁。 FSO742。未來科學網（Atypon），<https://doi.org/10.2144/fsoa-2021-0062>.                                                   |
| 激酶的親電彈頭               | 24  | Petri、László 等人。 “用於繪製蛋白質激酶中易處理半胱氨酸的反應性和可及性的親電彈頭庫。”歐洲藥物化學雜誌，第一卷。 207，2020 年 12 月，第112836.考研，<https://doi.org/10.1016/j.ejmech.2020.112836>.                                       |
| 激酶的特權支架               | 29  | 胡華斌等。 “競爭性和變構激酶抑製劑的系統比較揭示了共同的結構特徵。”歐洲藥物化學雜誌，第一卷。 214，2021 年 3 月，第113206. 科學直接,<https://doi.org/10.1016/j.ejmech.2021.113206>.                                                     |
| BRaf 抑製劑              | 54  | Agianian、Bogos 和 Evripidis Gavathiotis。 “BRAF 抑製劑在癌症中的最新見解。”藥物化學雜誌，第一卷。 61，沒有。 14，2018 年 7 月，第 5775-93 頁。 ACS出版物，<https://doi.org/10.1021/acs.jmedchem.7b01306>.                  |
| 常見氨基酸保護基團             | 346 | Isidro-Llobet，阿爾伯特等人。 “氨基酸保護基團。”化學評論，第一卷。 109，沒有。 6，2009 年 6 月，第 2455-504 頁。 DOI.org（交叉引用），<https://doi.org/10.1021/cr800323s>.                                                   |
| 新興的全氟烷基               | 27  | 佩爾奇，凱瑟琳 E.，等人。 “PFAS 健康影響數據庫：系統證據圖協議。”環境國際，第一卷。 130，2019 年 9 月，第104851. 科學直接,<https://doi.org/10.1016/j.envint.2019.05.045>.                                                      |
| 粘土吸附用化學品              | 33  | 奧爾，明日香 A.，等人。 “結合實驗等溫線、簡約模擬和模型來理解和預測蒙脫石粘土上的化學吸附。” ACS 歐米茄，第一卷。 6，沒有。 22，2021 年 6 月，第 14090-103 頁。考研，<https://doi.org/10.1021/acsomega.1c00481>.                                   |
| 大麻素                   | 63  | 特納，卡爾頓 E.，等人。 “大麻苜蓿 L. XVII 的成分。天然成分回顧。”天然產物雜誌，第一卷。 43，沒有。 2，1980 年 3 月，第 169-234 頁。 ACS出版物，<https://doi.org/10.1021/np50008a001>.                                                |
| 附表 1 美國麻醉品            | 240 | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                     |
| 附表 2 美國麻醉品            | 60  | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                     |
| 附表 3 美國麻醉品            | 22  | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                     |
| 附表 4 美國麻醉品            | 77  | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                     |
| 附表 5 美國麻醉品            | 8   | ECFR:: 21 CFR Part 1308 - 附表。                                                                                                                                                     |
| 皮卡爾                   | 179 | 舒爾金、亞歷山大 T. 和安舒爾金。 Pihkal：化學愛情故事。 1. ed.，8. 印刷，Transform，2010。                                                                                                                    |
| 輔料西咪替丁和阿昔洛韋           | 14  | Vaithianathan、Soundarya 等人。 “常見賦形劑對生物藥劑學分類系統第 3 類藥物西咪替丁和阿昔洛韋口服藥物吸收的影響。”藥學雜誌，第一卷。 105，沒有。 2，2016 年 2 月，第 996-1005 頁。考研，<https://doi.org/10.1002/jps.24643>.                        |
| 鎳雙膦酸鹽配體               | 不適用 | Clevenger、Andrew L. 等人。 “在鎳催化中使用二齒膦作為配體的趨勢。”化學評論，第一卷。 120，沒有。 13，2020 年 7 月，第 6124-96 頁。 DOI.org（交叉引用），<https://doi.org/10.1021/acs.chemrev.9b00682>.                             |
| 常見的正則表達式模式            | 1   |                                                                                                                                                                                   |

</details>                            

# GlobalChemExtensions

您可以使用多種工具來瀏覽和分析數據，並且可以在 google colab 演示或 Gitbook 文檔中找到不同應用程序的完整列表。使用 plotly 和 bokeh 設計的數據可視化擴展的演示如下所示：

<p align="center">
  <img width="800" height="600" src="https://raw.githubusercontent.com/Sulstice/global-chem/master/images/figures/figure_10.png">
</p>

<details><summary><h3>Extension List</h1><br/></summary>

| 延期              | 描述                                                                            | 應用                     |
| --------------- | ----------------------------------------------------------------------------- | ---------------------- |
| GlobalChem 化學實體 | GlobalChem 具有內部分子對象，具有關聯的所有常見屬性並轉換為 SMILES                                    | 力場                     |
| GlobalChem 生物實體 | GlobalChem 具有內部 DNA/RNA/蛋白質/分子對象，具有關聯的所有常見屬性並轉換為 SMILES                       | 生物信息學                  |
| 可視化 DNA/RNA 鏈   | 可視化 DNA 和 RNA 鏈並為其添加標籤                                                        | 生物信息學                  |
| 力場分子            | GlobalChem 可以將 CGenFF 和 GaFF2 文件作為對象解析、操作和寫入                                  | 力場                     |
| PDF 生成和解析       | GlobalChem 可以將 SMILES 生成 PDF 並將 PDF 轉換為 SMILES                                | 化學信息學                  |
| 微笑驗證            | GlobalChem 與 PySMILES、DeepSMILES、PartialSmiles、SELFIES、MolVS 連接，用於驗證 SMILES 集 | 化學信息學                  |
| 微笑質子狀態          | GlobalChem 可以採用一組化合物並預測 SMILES 串在一定 pH 值範圍內的質子化狀態                             | 化學信息學                  |
| 開源數據庫監控         | GlobalChem 使用 Uptime-Cheminformatics 跟踪開源化學數據                                 | development_operations |
| Networkx 軟件適配器  | GlobalChem Network 可以轉換為 NetworkX Graph Objects                               | 化學信息學                  |
| SMARTS 模式驗證     | GlobalChem 使用 MiniFrag 數據庫測試 SMARTS 字符串在功能組選擇中的準確性                            | 化學信息學                  |
| 主成分分析           | GlobalChem 可以輕鬆解釋 SMILES、指紋、聚類和應用 PCA 分析，用戶可以調整參數                             | 化學信息學                  |
| 藥物設計過濾器         | GlobalChem 可以根據通用藥物設計過濾規則過濾化合物                                                | 化學信息學                  |
| 深層散射分析          | 為了可視化分子組之間的關係，GlobalChem 提供了平行坐標圖生成                                           | 化學信息學                  |
| 旭日形徑向分析         | GlobalChem 提供了一種旭日機制，允許用戶觀察化合物組與公共組的關係                                        | 化學信息學                  |
| 圖形模板            | GlobalChem 提供圖形模板以幫助更快地進行數據分析，目前唯一提供的是 Plotly                                 | 化學信息學                  |
| CGenFF 差異分數     | GlobalChem 可以根據原子類型提供兩種分子之間的差異                                                | 力場                     |
| 一種熱編碼           | GlobalChem 擁有自己的一個基於機器學習通用列表的熱門編碼器和解碼器                                        | 化學信息學                  |
| SMARTS 模式標識符    | GlobalChem 連接到 SMARTS Plus，可以提供不同 SMARTS 組件的可視化                               | 化學信息學                  |
| Psi4 解析器        | 提供 Psi4 輸出文件的解析和提取值                                                           | 量子化學                   |
| 坐標商店            | 用於在 xyz 和 z 矩陣中分佈的小分子坐標的倉庫                                                    | 量子化學                   |
| 可視化分子軌道         | 可視化來自 Psi4 輸出 cubeprop 的多維數據集文件                                               | 量子化學                   |

</details>

# 開源軟件合規性

`GlobalChem`遵循聯邦法規第 21 篇第 11 部分概述的相同原則；電子記錄，
電子簽名（21 CFR 第 11 部分）指導文件。由於沒有關於如何處理開源軟件的正式指南，我們
嘗試完成要求。 FDA 認為第 11 部分適用於以下電子記錄標準以及如何`GlobalChem`完成每個組件：

-   **合理性：**`GlobalChem`建立在使用閱讀和重繪從書籍和論文中提取的數據之上。它添加了一個組件
    IUPAC/SMILES/SMARTS 字符串以電子方式存儲它，從而為它的數據提供獨特的組件。記錄是開源的
    並由存儲庫的維護者和開源社區反饋進行適當的版本控制。`GlobalChem`的目的仍然未知，因為它進入了開源部署。我們已經構建了擴展功能
    一個單獨的包`GlobalChemExtensions`這確實取決於`GlobalChem`.由於每個版本都被適當地打包，
    需要依賴一個版本，然後它的軟件才可用`Github`和`PyPi`.標準操作程序 (SOP)
    可以從維護的擴展實用程序文檔中提交`Gitbook`.

-   **驗證：**`GlobalChem`遵循良好自動化製造規範 (GAMP) 類別 3，即“安裝後使用的軟件”
    並且可能是“可配置的”。`GlobalChem`測試來自內部，文檔作為最終測試
    功能，因為這是用戶測試最多的東西，因為我們依賴開源。持續集成 (CI)
    系統還同時構建以作為基本功能測試`GlobalChem`圖網絡。存儲的數據
    由該領域的專家維護，但如果發現錯誤，可能會根據社區反饋進行更改。

-   **審計追踪：**`GlobalChem`受版本控制`Git`並託管在微軟的平台上`Github`.`GlobalChem`遵循語義
    模式的版本控制`X1.X2.X3`:`X1`用測試、文檔和平均值標記正式的穩定版本
    對軟件或功能進行大的重構，`X2`意味著添加了一個新功能，無論是否有測試和文檔，但
    如此迭代。`X3`表示要添加到函數中的“熱”修復（一個簡單的錯誤）、小功能或附加參數
    ，或迭代數據。

-   **舊系統：**`GlobalChem`自第一個版本發布以來已經運行了近 2 年`0.3.0`2020 年 5 月。`GlobalChem`在開源社區中建立了完整的路徑，每個版本都被編目並且對所有人可見。這滿足
    規則概述了確定遺留系統。我們使用社交媒體平台（Twitter、Github、LinkedIn）提供的社區反饋
    作為書面證據和理由`GlobalChem`適合化學信息學的預期用途。

-   **記錄副本：**`GlobalChem`有記錄存儲在`Github`對於可以導出為各種格式的軟件
    微軟。對於文檔，它託管在`Gitbook`並根據軟件進行版本控制。每一個“書”
    可以導出為適用於 FDA 提交的便攜式數據格式 (PDF)。

-   **記錄保留：**`GlobalChem`具有控制為唯一 ID (UUID) 的文檔版本記錄，該唯一 ID (UUID) 用作其標識符
    對於存儲在`Gitbook`.每個版本都存儲為降價文件，並在需要時轉換為 PDF。

`GlobalChem`具有 Mozilla 公共許可證 2.0 版。`GlobalChem`允許您在更大的工作中使用該軟件，並且
如果您願意，可以通過修改對其進行擴展。意外情況是，如果您安裝`GlobalChem`並發布新軟件
那麼您必須遵循我們的開源社區許可證中安裝的相同原則。

# 數據採集

參考文獻和相關化合物列表是根據科學貢獻者的興趣選擇的。這應包括考慮與科學界的相關性。
SMILES 字符串可以用多種方法進行抽象：

-   對於簡單的分子，可以使用視覺直接翻譯 SMILES 的一種表示
    檢查。這通常適用於報告列表開頭包含最常見分母環的化合物。

-   對於復雜的分子，圖像可以在免費版本的 ChemDraw 中重新繪製，然後翻譯成 SMILES。

-   對於寫入 SMILES 且不知道 IUPAC 的來源，SMILES 將被翻譯成 ChemDraw 並檢索名稱。
    請注意，某些名稱可能會根據人工檢查進行修改，以支持首選名稱。

-   對於聚合物紙，名稱中省略了站點點，並且根據首選名稱調整了一些命名法
    超過傳統。例如：用於標記聚合物連接點的“yl”已被刪除，以降低英語的複雜性。

-   在自由基的情況下，一些 SMILES 被調整以去除自由基化學特徵，因為它們用作連接點。然而，在某些情況下，自由基成分得以保留，特別是在 IUPAC 藍皮書常見取代基的情況下。

-   SMARTS 字符串是使用 RDKit (4) 從 SMILES 改編而來的。

* * *

# 引文

它正在路上

# 許可

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FSulstice%2Fglobal-chem?ref=badge_large)
