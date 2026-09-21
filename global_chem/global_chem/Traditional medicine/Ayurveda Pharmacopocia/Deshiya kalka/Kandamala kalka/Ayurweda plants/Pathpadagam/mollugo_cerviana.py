#!/usr/bin/env python3
#
# GlobalChem - Mollugo Cerviana
#
# -----------------------------

# Mollugo cerviana (L.) Ser. - whole plant chemical library
#
# Ayurvedic name  : Pathpadagam
# English names   : Threadstem, Slender Carpetweed, Carpetweed
# Family          : Molluginaceae
# Plant material  : Whole plant (30 compounds)
#
# Evidence tiers
# --------------
# T1 : reported for Mollugo cerviana itself (refs 1-3)
# T2 : GC-MS peak of a combined ethanolic extract of Mollugo cerviana and
#      Mukia maderaspatana (ref 4)
#
# | # | Common Name | Class | Formula | Tier |
# |---|-------------|-------|---------|------|
# | 1 | orientin | flavone C-glycoside | C21H20O11 | T1 |
# | 2 | vitexin | flavone C-glycoside | C21H20O10 | T1 |
# | 3 | quercetin | flavonol | C15H10O7 | T1 |
# | 4 | gallic acid | phenolic acid | C7H6O5 | T1 |
# | 5 | caffeic acid | hydroxycinnamic acid | C9H8O4 | T1 |
# | 6 | stigmasterol | phytosterol | C29H48O | T1 |
# | 7 | beta-sitosterol | phytosterol | C29H50O | T1 |
# | 8 | lupeol | lupane triterpenoid | C30H50O | T1 |
# | 9 | ursolic acid | ursane triterpenoid | C30H48O3 | T1 |
# | 10 | ascorbic acid | vitamin C | C6H8O6 | T1 |
# | 11 | alpha-tocopherol | vitamin E | C29H50O2 | T1 |
# | 12 | retinol | vitamin A | C20H30O | T1 |
# | 13 | colchicine | tropolone alkaloid | C22H25NO6 | T1 |
# | 14 | ergost-7-en-3-ol | phytosterol | C28H48O | T2 |
# | 15 | phytol | diterpene alcohol | C20H40O | T2 |
# | 16 | neophytadiene | diterpene hydrocarbon | C20H38 | T2 |
# | 17 | 6,10,14-trimethylpentadecan-2-one | isoprenoid ketone | C18H36O | T2 |
# | 18 | ethyl palmitate | fatty acid ester | C18H36O2 | T2 |
# | 19 | alpha-linolenic acid | fatty acid | C18H30O2 | T2 |
# | 20 | ethyl linoleate | fatty acid ester | C20H36O2 | T2 |
# | 21 | ethyl linolenate | fatty acid ester | C20H34O2 | T2 |
# | 22 | methyl 17-methyloctadecanoate | fatty acid ester | C20H40O2 | T2 |
# | 23 | 2-palmitoylglycerol | monoacylglycerol | C19H38O4 | T2 |
# | 24 | glyceryl monostearate | monoacylglycerol | C21H42O4 | T2 |
# | 25 | squalene | triterpene hydrocarbon | C30H50 | T2 |
# | 26 | gamma-tocopheryl methyl ether | vitamin E derivative | C29H50O2 | T2 |
# | 27 | ascorbyl 2,6-dipalmitate | vitamin C ester | C38H68O8 | T2 |
# | 28 | 2,5-di-tert-butylhydroquinone | hydroquinone | C14H22O2 | T2 |
# | 29 | diethyl phthalate | phthalate ester | C12H14O4 | T2 |
# | 30 | di-n-octyl phthalate | phthalate ester | C24H38O4 | T2 |
#
# References
# ----------
# 1. Padmapriya, P., and S. Maneemegalai. "Qualitative and Quantitative Analysis
#    of the Phytochemical Constituents of Mollugo cerviana (L.)." International
#    Journal of Pharmacy and Drug Analysis, vol. 2, no. 9, 2014, pp. 695-699.
# 2. Arthi, P., et al. "An Integrative Review of Mollugo cerviana: From Weed to
#    Medicinal Plant." International Journal of Progressive Research in
#    Engineering Management and Science, vol. 5, no. 8, 2025, pp. 573-578.
# 3. Rebecca, F. S., R. P. Nivedha, and R. Sharmila. "Molecular Docking Analysis
#    of Bioactive Compounds from Mollugo cerviana (L.) Ser with DHFR for
#    Antifungal Activity." Bioinformation, vol. 17, no. 11, 2021, pp. 944-948.
#    https://doi.org/10.6026/97320630017944
# 4. Rebecca, F. S., et al. "Antimicrobial Evaluation and Molecular Docking
#    Studies of the Combined Ethanolic Extract of Mollugo cerviana and Mukia
#    maderaspatana." Bioinformation, vol. 19, no. 2, 2023, pp. 190-195.
#    https://doi.org/10.6026/97320630019190
#
# Contributor: Janith Madhusanka (GitHub: janithmadhu277-cmyk)
# SMARTS strings were adapted from the SMILES using the RDKit atomic-number
# SMARTS convention and validated against the SMILES with RDKit.

class MollugoCerviana(object):

    def __init__(self):

        self.name = 'mollugo_cerviana'

    @staticmethod
    def get_smiles():

        smiles = {
            'orientin': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O',
            'vitexin': 'C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O',
            'quercetin': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O',
            'gallic acid': 'C1=C(C=C(C(=C1O)O)O)C(=O)O',
            'caffeic acid': 'C1=CC(=C(C=C1/C=C/C(=O)O)O)O',
            'stigmasterol': 'CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C',
            'beta-sitosterol': 'CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C',
            'lupeol': 'CC(=C)[C@@H]1CC[C@@]2([C@H]1[C@H]3CC[C@@H]4[C@@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C',
            'ursolic acid': 'C[C@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H](C5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C(=O)O',
            'ascorbic acid': 'C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O',
            'alpha-tocopherol': 'CC1=C(C2=C(CC[C@@](O2)(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=C1O)C)C',
            'retinol': 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C',
            'colchicine': 'CC(=O)N[C@H]1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC',
            'ergost-7-en-3-ol': 'CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2=CCC4C3(CCC(C4)O)C)C',
            'phytol': 'CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC/C(=C/CO)/C',
            'neophytadiene': 'CC(C)CCCC(C)CCCC(C)CCCC(=C)C=C',
            '6,10,14-trimethylpentadecan-2-one': 'CC(=O)CCCC(C)CCCC(C)CCCC(C)C',
            'ethyl palmitate': 'CCCCCCCCCCCCCCCC(=O)OCC',
            'alpha-linolenic acid': 'CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O',
            'ethyl linoleate': 'CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)OCC',
            'ethyl linolenate': 'CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)OCC',
            'methyl 17-methyloctadecanoate': 'COC(=O)CCCCCCCCCCCCCCCC(C)C',
            '2-palmitoylglycerol': 'CCCCCCCCCCCCCCCC(=O)OC(CO)CO',
            'glyceryl monostearate': 'CCCCCCCCCCCCCCCCCC(=O)OCC(O)CO',
            'squalene': 'CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC=C(C)C',
            'gamma-tocopheryl methyl ether': 'COC1=C(C)C(C)=C2OC(C)(CCCC(C)CCCC(C)CCCC(C)C)CCC2=C1',
            'ascorbyl 2,6-dipalmitate': 'CCCCCCCCCCCCCCCC(=O)OC[C@@H]([C@@H]1C(=C(C(=O)O1)OC(=O)CCCCCCCCCCCCCCC)O)O',
            '2,5-di-tert-butylhydroquinone': 'CC(C)(C)C1=CC(=C(C=C1O)C(C)(C)C)O',
            'diethyl phthalate': 'CCOC(=O)C1=CC=CC=C1C(=O)OCC',
            'di-n-octyl phthalate': 'CCCCCCCCOC(=O)C1=CC=CC=C1C(=O)OCCCCCCCC',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'orientin': '[#8]-[#6]-[#6]1-[#8]-[#6](-[#6]2:[#6](-[#8]):[#6]:[#6](-[#8]):[#6]3:[#6]:2:[#8]:[#6](-[#6]2:[#6]:[#6]:[#6](-[#8]):[#6](-[#8]):[#6]:2):[#6]:[#6]:3=[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]',
            'vitexin': '[#8]-[#6]-[#6]1-[#8]-[#6](-[#6]2:[#6](-[#8]):[#6]:[#6](-[#8]):[#6]3:[#6]:2:[#8]:[#6](-[#6]2:[#6]:[#6]:[#6](-[#8]):[#6]:[#6]:2):[#6]:[#6]:3=[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]',
            'quercetin': '[#8]=[#6]1:[#6](-[#8]):[#6](-[#6]2:[#6]:[#6]:[#6](-[#8]):[#6](-[#8]):[#6]:2):[#8]:[#6]2:[#6]:[#6](-[#8]):[#6]:[#6](-[#8]):[#6]:1:2',
            'gallic acid': '[#8]=[#6](-[#8])-[#6]1:[#6]:[#6](-[#8]):[#6](-[#8]):[#6](-[#8]):[#6]:1',
            'caffeic acid': '[#8]=[#6](-[#8])-[#6]=[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]):[#6](-[#8]):[#6]:1',
            'stigmasterol': '[#6]-[#6]-[#6](-[#6]=[#6]-[#6](-[#6])-[#6]1-[#6]-[#6]-[#6]2-[#6]-1(-[#6]-[#6]-[#6]3-[#6]-2-[#6]-[#6]=[#6]4-[#6]-3(-[#6]-[#6]-[#6](-[#6]-4)-[#8])-[#6])-[#6])-[#6](-[#6])-[#6]',
            'beta-sitosterol': '[#6]-[#6]-[#6](-[#6]-[#6]-[#6](-[#6])-[#6]1-[#6]-[#6]-[#6]2-[#6]-1(-[#6]-[#6]-[#6]3-[#6]-2-[#6]-[#6]=[#6]4-[#6]-3(-[#6]-[#6]-[#6](-[#6]-4)-[#8])-[#6])-[#6])-[#6](-[#6])-[#6]',
            'lupeol': '[#6]-[#6](=[#6])-[#6]1-[#6]-[#6]-[#6]2(-[#6]-1-[#6]3-[#6]-[#6]-[#6]4-[#6]5(-[#6]-[#6]-[#6](-[#6](-[#6]-5-[#6]-[#6]-[#6]-4(-[#6]-3(-[#6]-[#6]-2)-[#6])-[#6])(-[#6])-[#6])-[#8])-[#6])-[#6]',
            'ursolic acid': '[#6]-[#6]1-[#6]-[#6]-[#6]2(-[#6]-[#6]-[#6]3(-[#6](=[#6]-[#6]-[#6]4-[#6]-3(-[#6]-[#6]-[#6]5-[#6]-4(-[#6]-[#6]-[#6](-[#6]-5(-[#6])-[#6])-[#8])-[#6])-[#6])-[#6]-2-[#6]-1-[#6])-[#6])-[#6](=[#8])-[#8]',
            'ascorbic acid': '[#6](-[#6](-[#6]1-[#6](=[#6](-[#6](=[#8])-[#8]-1)-[#8])-[#8])-[#8])-[#8]',
            'alpha-tocopherol': '[#6]-[#6]1(-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6])-[#6]-[#6]-[#6]2:[#6](-[#6]):[#6](-[#8]):[#6](-[#6]):[#6](-[#6]):[#6]:2-[#8]-1',
            'retinol': '[#6]-[#6]1=[#6](-[#6](-[#6]-[#6]-[#6]-1)(-[#6])-[#6])-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6](=[#6]-[#6]-[#8])-[#6])-[#6]',
            'colchicine': '[#6]-[#8]-[#6]1:[#6]:[#6]2:[#6](:[#6](-[#8]-[#6]):[#6]:1-[#8]-[#6])-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6](=[#8]):[#6]:[#6]:1-[#6](-[#7]-[#6](-[#6])=[#8])-[#6]-[#6]-2',
            'ergost-7-en-3-ol': '[#6]-[#6](-[#6])-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]1-[#6]-[#6]-[#6]2-[#6]-1(-[#6]-[#6]-[#6]3-[#6]-2=[#6]-[#6]-[#6]4-[#6]-3(-[#6]-[#6]-[#6](-[#6]-4)-[#8])-[#6])-[#6]',
            'phytol': '[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](=[#6]-[#6]-[#8])-[#6]',
            'neophytadiene': '[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](=[#6])-[#6]=[#6]',
            '6,10,14-trimethylpentadecan-2-one': '[#6]-[#6](=[#8])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            'ethyl palmitate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6]',
            'alpha-linolenic acid': '[#6]-[#6]-[#6]=[#6]-[#6]-[#6]=[#6]-[#6]-[#6]=[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]',
            'ethyl linoleate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]=[#6]-[#6]-[#6]=[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6]',
            'ethyl linolenate': '[#6]-[#6]-[#6]=[#6]-[#6]-[#6]=[#6]-[#6]-[#6]=[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6]',
            'methyl 17-methyloctadecanoate': '[#6]-[#8]-[#6](=[#8])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '2-palmitoylglycerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6](-[#6]-[#8])-[#6]-[#8]',
            'glyceryl monostearate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6](-[#8])-[#6]-[#8]',
            'squalene': '[#6]-[#6](-[#6])=[#6]-[#6]-[#6]-[#6](-[#6])=[#6]-[#6]-[#6]-[#6](-[#6])=[#6]-[#6]-[#6]-[#6]=[#6](-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6]',
            'gamma-tocopheryl methyl ether': '[#6]-[#6]1(-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6])-[#6]-[#6]-[#6]2:[#6]:[#6](-[#8]-[#6]):[#6](-[#6]):[#6](-[#6]):[#6]:2-[#8]-1',
            'ascorbyl 2,6-dipalmitate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6](-[#6]1-[#6](=[#6](-[#6](=[#8])-[#8]-1)-[#8]-[#6](=[#8])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])-[#8])-[#8]',
            '2,5-di-tert-butylhydroquinone': '[#6]-[#6](-[#6])(-[#6])-[#6]1:[#6]:[#6](-[#8]):[#6](-[#6](-[#6])(-[#6])-[#6]):[#6]:[#6]:1-[#8]',
            'diethyl phthalate': '[#6]-[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](=[#8])-[#8]-[#6]-[#6]',
            'di-n-octyl phthalate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](=[#8])-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
        }

        return smarts
