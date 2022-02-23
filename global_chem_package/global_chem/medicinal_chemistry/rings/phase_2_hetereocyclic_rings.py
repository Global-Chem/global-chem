#!/usr/bin/env python3
#
# GlobalChem - Phase 2 Hetereocyclic Rings
#
# ----------------------------------------

class Phase2HetereoCyclicRings(object):

    def __init__(self):

        self.name = 'phase2_hetereo_cyclic_rings'

    @staticmethod
    def get_smiles():

        smiles = {
            'pyridine': 'C1=CC=NC=C1',
            'indole': 'C12=CC=CC=C1C=CN2',
            'imidazole': 'C1=CN=CN1',
            'thiazol-2-amine': 'NC1=NC=CS1',
            'tetrazole': 'C1=NN=NN1',
            '1,2,4-triazole': 'C1=NC=NN1',
            'thiophene': 'C1=CC=CS1',
            'cytosine': 'O=C1N=C(N)C=CN1',
            'adenine': 'NC1=NC=NC2=C1N=CN2',
            '5-methylindole': 'CC1=CC=C2C(C=CN2)=C1',
            'isocaffeine': 'O=C(N1C)NC2=C(N=CN2)C1=O',
            'tetrazolethiol': 'SN1N=NN=C1',
            '3-methylisoxazole': 'C1=CC=NO1',
            '1-methylimidazole': 'CN1C=NC=C1',
            '2-methylimidazole': 'CC1=NC=CN1',
            'guanine': 'NC(N1)=NC2=C(N=CN2)C1=O',
            'quinoline': 'C12=CC=CC=C1N=CC=C2',
            'furan': 'C1=CC=CO1',
            'tosufloxacin': 'NC1=C(F)C=C2C(NC=C(C(O)=O)C2=O)=N1'
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'pyridine': '[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1',
            'indole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#7H]:2',
            'imidazole': '[#6]1:[#6]:[#7]:[#6]:[#7H]:1',
            'thiazol-2-amine': '[#7]-[#6]1:[#7]:[#6]:[#6]:[#16]:1',
            'tetrazole': '[#6]1:[#7]:[#7]:[#7]:[#7H]:1',
            '1,2,4-triazole': '[#6]1:[#7]:[#6]:[#7]:[#7H]:1',
            'thiophene': '[#6]1:[#6]:[#6]:[#6]:[#16]:1',
            'cytosine': '[#8]=[#6]1:[#7]:[#6](-[#7]):[#6]:[#6]:[#7H]:1',
            'adenine': '[#7]-[#6]1:[#7]:[#6]:[#7]:[#6]2:[#6]:1:[#7]:[#6]:[#7H]:2',
            '5-methylindole': '[#6]-[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:[#6]:[#7H]:2):[#6]:1',
            'isocaffeine': '[#8]=[#6]1:[#7](-[#6]):[#6](:[#6]2:[#6](:[#7H]:1):[#7H]:[#6]:[#7]:2)=[#8]',
            'tetrazolethiol': '[#16]-[#7]1:[#7]:[#7]:[#7]:[#6]:1',
            '3-methylisoxazole': '[#6]1:[#6]:[#6]:[#7]:[#8]:1',
            '1-methylimidazole': '[#6]-[#7]1:[#6]:[#7]:[#6]:[#6]:1',
            '2-methylimidazole': '[#6]-[#6]1:[#7]:[#6]:[#6]:[#7H]:1',
            'guanine': '[#7]-[#6]1:[#7H]:[#6](:[#6]2:[#6](:[#7]:1):[#7H]:[#6]:[#7]:2)=[#8]',
            'quinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6]:[#6]:[#6]:2',
            'furan': '[#6]1:[#6]:[#6]:[#6]:[#8]:1',
            'tosufloxacin': '[#7]-[#6]1:[#6](-[#9]):[#6]:[#6]2:[#6](:[#7H]:[#6]:[#6](-[#6](-[#8])=[#8]):[#6]:2=[#8]):[#7]:1',
        }

        return smarts
