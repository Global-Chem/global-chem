#!/usr/bin/env python3
#
# GlobalChem - Insect Sex Pheromones
#
# ------------------------------------

class InsectSexPheromones(object):

    def __init__(self):

        self.name = 'insect_sex_pheromones'

    @staticmethod
    def get_smiles():
        smiles = {
            '2,2-dimethyl-3-isopropylidene cyclopropyl': 'C1(C)(C)=C(OC(=O)CC)C1(C)(C)',
            'cis-7-dodecen-1-ol acetate': 'CCCC/C=C\CCCCCCOC(=O)C',
            'cis-9-tetradecen-1-ol acetate': 'CCCC/C=C\CCCCCCCCOC(=O)C',
            'trans-2-hexen-1-ol acetate': 'CCC/C=C/COC(=O)C',
            'trans-2-hexen-1-ol butyrate': 'CCC/C=C/COC(=O)CCC',
            'trans-2-octen-1-ol acetate': 'CCCCC/C=C/COC(=O)C',
            'trans-2-decen-1-ol acetate': 'CCCCCCC/C=C/COC(=O)C',
            'undecanal': 'CCCCCCCCCCC=O',
            'cis-11-octadecenal': 'CCCCCC\C=C/CCCCCCCCCC=O',
            '11-acetoxyundecanal': 'O=CCCCCCCCCCCOC(=O)C',
            '9-acetoxynonanal': 'O=CCCCCCCCCOC(=O)C',
            'cis-9,trans-12-tetradecadien-1-ol acetate': 'C/C=C/C/C=C\CCCCCCCCOC(=O)C',
            'cis-11-tetradecen-1-ol acetate': 'CC\C=C/CCCCCCCCCCOC(=O)C',
            'trans-10,cis-12-hexadecadien-1-ol': 'CCC\C=C/C=C/CCCCCCCCCO',
            'trans-11-tetradecenal': 'O=CCCCCCCCCC/C=C/CC',
            'cis-11-tetradecen-1-ol acetate': 'CC\C=C/CCCCCCCCCCOC(=O)C',
            'trans-7-dodecen-1-ol acetate': 'CCCC/C=C/CCCCCCOC(=O)C',
            'trans-2,trans-6-decadien-1,10-diol': 'OC/C=C/CC/C=C/CCCO',
            '10-hydroxy-3,7-dimethyl-trans-2,trans-6 decadienoic acid': 'OCCCC/C(C)=C/CC/C(C)=C/C(=O)O',
            'cis-8-dodecen-1-ol acetate': 'CCC/C=C\CCCCCCCOC(=O)C',
            '2-methylheptadecane': 'CC(C)CCCCCCCCCCCCCCC',
            'trans-8,trans-10-dodecadien-1-ol': 'OCCCCCCC/C=C/C=C/C',
            '2,3-dihydro-7-methyl-1H-pyrrolizin-1-one': 'O=C1CCN2C=CC(C)=C12',
            'n-hexadecyl acetate': 'CCCCCCCCCCCCCCCCOC(=O)C',
            'cis-11-octadecen-1-ol acetate': 'CCCCCC\C=C/CCCCCCCCCCOC(=O)C',
            '10-propyl-trans-5,9-tridecadien-1-ol acetate': 'CCC/C(CCC)=C/CC/C=C/CCCCOC(=O)C',
            'N,N-diethyl-m-toluamide': 'N(CC)(CC)C(=O)C1=CC(C)=CC=C1',
            'd-10-acetoxy-cis-7-hexadecen-1-ol': 'OCCCCCC\C=C/CC(OC(=O)C)CCCCCC',
            'cis-7,8-epoxy-2-methyloctadecone': 'CCCCCCCCCC[C@H]1O[C@H]1CCCCC(C)C',
            'cis-7-dodecen-1-ol acetate': 'CCCC/C=C\CCCCCCOC(=O)C',
            '(-)-methyl trans-2,4,5-tetradecatrienoate': 'CCCCCCCCC=C=C/C=C/C(=O)OC',
            'trans-3,cis-5-tetradecadienoic acid': 'CCCCCCCCC=CC=CCC(=O)O',
            'exo-7-ethyl-5-methyl-6,8-dioxabicyclo[3.2.1]octane': 'O1[C@H](CC)C(O2)CCCC12C',
            '1,5-dimethyl-6.8-dioxabicyclo[3.2.1]octane': 'O1CC(C)(O2)CCCC12C',
            'trans-verbenol': 'C1(C)(C)C2CC1C(O)C=C2C',
            '(-)-14-methyl-cis-8-hexadecen-1-ol': 'CC[C@@H](C)CCCC/C=C\CCCCCCCO',
            '(-)-methyl 14-methyl-cis-8-hexadecenoate': 'CCC(C)CCCC/C=C\CCCCCCC(=O)OC'
        }

        return smiles