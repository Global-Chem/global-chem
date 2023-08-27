#!/usr/bin/env python3
#
# GlobalChem - Gin
#
# -----------------------------------

class Gin(object):

    def __init__(self):

        self.name = 'gin'

    @staticmethod
    def get_smiles():

        '''

        Missing:
        'verbeneen': 
        'I-phellandrene'
        "t-ocimene"
        "c-rose oxide"
        "c-sabinene hydrate"
        "l-4-terpine"
        "terpenyl acetate"
        "t-beta-caryophyllene"
        "sesquiterpene ni"
        "t-beta-farnesene"
        "t-muurulol" OR "T-Muurolol"
        '''

     

        smiles = {
            'alpha-pinene': "CC1=CCC2CC1C2(C)C",
            'alpha-thujene': "CC1=CCC2(C1C2)C(C)C",
            'camphene': "CC1(C2CCC(C2)C1=C)C",
            'beta-pinene': "CC1(C2CCC(=C)C1C2)C",
            "sabinene": "CC(C)C12CCC(=C)C1C2",
            'delta-3-carene': 'CC1=CCC2C(C1)C2(C)C',
            "beta-Myrcene": "CC(=CCCC(=C)C=C)C",
            "alpha-terpinene": "CC1=CC=C(CC1)C(C)C",
            "DL-limonene": "CC1=CCC(CC1)C(=C)C",
            "beta-phellandrene": "CC(C)C1CCC(=C)C=C1",
            "gamma-terpinene": "CC1=CCC(=CC1)C(C)C",
            "p-cymene": "CC1=CC=C(C=C1)C(C)C",
            "alpha-terpinolene": "CC1=CCC(=C(C)C)CC1",
            "verbenyl ethyl ether": "CCOC1C=C(C2CC1C2(C)C)C",
            "citronellal": "CC(CCC=C(C)C)CC=O",
            "campholenal": "CC1=CC[C@H](C1(C)C)CC=O",
            "camphor": "CC1(C2CCC1(C(=O)C2)C)C",
            "linalool": "CC(=CCCC(C)(C=C)O)C",
            "bornyl acetate": "CC(=O)O[C@@H]1C[C@@H]2CC[C@]1(C2(C)C)C",
            "myrtenal": "CC1(C2CC=C(C1C2)C=O)C",
            "alpha-terpineol": "CC1=CCC(CC1)C(C)(C)O",
            "neryl acetate": "CC(=CCC/C(=C\COC(=O)C)/C)C",
            "geranyl acetate": "CC(=CCC/C(=C/COC(=O)C)/C)C",
            "cuminal": "CC(C)C1=CC=C(C=C1)C=O",
            "beta-citronellol": "CC(CCC=C(C)C)CCO",
            "myrtenol": "CC1(C2CC=C(C1C2)CO)C",
            "nerol": "CC(=CCC/C(=C\CO)/C)C",
            "t-carveol": "CC1=CC[C@H](C[C@@H]1O)C(=C)C",
            "t-geraniol": "CC(=CCC/C(=C/CO)/C)C",
            "alpha-cubebene": "C[C@@H]1CC[C@H]([C@H]2[C@]13[C@@H]2C(=CC3)C)C(C)C",
            "alpha-copaene": "CC1=CCC2C3C1C2(CCC3C(C)C)C",
            "beta-cubebene": "C[C@@H]1CC[C@H]([C@H]2[C@]13[C@@H]2C(=C)CC3)C(C)C",
            "beta-elemene": "CC(=C)[C@@H]1CC[C@@]([C@@H](C1)C(=C)C)(C)C=C",
            "gamma-elemene": "CC(=C1CC[C@@]([C@@H](C1)C(=C)C)(C)C=C)C",
            "alpha-humulene": "C/C/1=C\CC(/C=C/C/C(=C/CC1)/C)(C)C",
            "gamma-muurolene": "CC1=C[C@@H]2[C@H](CC1)C(=C)CC[C@H]2C(C)C",
            "germacrene D": "C/C/1=C\CCC(=C)/C=C/[C@@H](CC1)C(C)C",
            "alpha-selinene": "CC1=CCC[C@]2([C@H]1C[C@@H](CC2)C(=C)C)C",
            "alpha-muurolene": "CC1=C[C@@H]2[C@H](CC1)C(=CC[C@H]2C(C)C)C",
            "delta-cadinene": "CC1=C[C@H]2[C@@H](CCC(=C2CC1)C)C(C)C",
            "gamma-cadinene": "CC1=C[C@H]2[C@H](CC1)C(=C)CC[C@@H]2C(C)C",
            "cadina-1,4-diene": "C[C@H]1CCC(C2C1=CCC(=C2)C)C(C)C",
            "germacrene B": "C/C/1=C\CC/C(=C/CC(=C(C)C)CC1)/C",
            "alpha-calacorene": "CC1=CC[C@H](C2=C1C=CC(=C2)C)C(C)C",
            "caryophyllene oxide": "C[C@@]12CC[C@@H]3[C@H](CC3(C)C)C(=C)CC[C@H]1O2",
            "torreyol": "CC1=C[C@H]2[C@@H](CC[C@@]([C@H]2CC1)(C)O)C(C)C",
            "elemol": "CC(=C)[C@@H]1C[C@@H](CC[C@@]1(C)C=C)C(C)(C)O",
            "spathulenol": "C[C@@]1(CC[C@@H]2[C@@H]1[C@H]3[C@H](C3(C)C)CCC2=C)O",
            "t-cadinol": "CC1=C[C@H]2[C@@H](CC[C@]([C@@H]2CC1)(C)O)C(C)C",
            "eudesmol": "C[C@@H]1CCC[C@]2([C@H]1C[C@@H](CC2)C(C)(C)O)C",
            "alpha-cadinol": "CC1=C[C@H]2[C@@H](CC[C@@]([C@@H]2CC1)(C)O)C(C)C",
            "nonanal": "CCCCCCCCC=O",
            "benzaldehyde": "C1=CC=C(C=C1)C=O",
            "2-undecanone": "CCCCCCCCCC(=O)C",
            "z-citral": "CC(=CCC/C(=C/C=O)/C)C",
            "Carvone": "CC1=CCC(CC1=O)C(=C)C",
            "Verbenone": "CC1=CC(=O)C2CC1C2(C)C",
            "trans-pinocarveol": "CC1([C@@H]2C[C@H]1C(=C)[C@H](C2)O)C",
            "hexanal": "CCCCCC=O",
            "heptanal": "CCCCCCC=O",
            "dehydrosabinene": "CC(C)C12CC1C(=C)C=C2",
            "limonene": "CC1=CCC(CC1)C(=C)C",
            "bicycloelemene": "CC(=C)[C@H]1[C@H]2[C@H](C2(C)C)CC[C@@]1(C)C=C",
        

        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            
        }

        return smarts