#!/usr/bin/env python3
#
# GlobalChem - PrivilegedScaffolds
#
# -----------------------------------

class PrivilegedScaffolds(object):

    def __init__(self):

        self.name = 'privileged_scaffolds'

    @staticmethod
    def get_smiles():

        smiles = {
            'indole': 'C12=CC=CC=C1C=CN2',
            'quinoline': 'C12=CC=CC=C1N=CC=C2',
            'isoquinoline': 'C12=C(C=NC=C2)C=CC=C1',
            'purine': 'C12=NC=NC=C1NC=N2',
            'quinoxaline': 'C12=CC=CC=C1N=CC=N2',
            'quinazolinone': 'O=C1NC2=C(C=CC=C2)C=N1',
            'tetrahydroisoquinoline': 'C12=C(CNCC2)C=CC=C1',
            'tetrahydraquinoline': 'C12=C(NCCC2)C=CC=C1',
            'benzoxazole': 'C12=CC=CC=C1OC=N2',
            'benzofuran': 'C12=CC=CC=C1C=CO2',
            '3,3-dimethylbenzopyran': 'CC1(C)C=CC2=CC=CC=C2O1',
            'chromone': 'O=C1C=COC2=C1C=CC=C2',
            'coumarin': 'O=C1OC2=C(C=CC=C2)C=C1',
            'carbohydrate': 'OCC1OC(O)C(O)C(O)C1O',
            'steroid': 'C12CCCCC1C3C(C(CCC4)C4CC3)CC2',
            'prostanoic acid': 'CCCCCCCC[C@@H]1[C@H](CCC1)CCCCCCC(O)=O',
            'benzodiazepine': 'O=C1CN=C(C2=CC=CC=C2)C3=C(C=CC=C3)N1',
            'arylpiperidine': 'C1(C2CCNCC2)=CC=CC=C1',
            'arylpiperizine': 'C1(N2CCNCC2)=CC=CC=C1',
            'benzylpiperidine': 'N1(CC2=CC=CC=C2)CCCCC1',
            'benzothiophene': 'C12=CC=CC=C1C=CS2',
            'dihydropyridine': 'C1CC=CC=N1',
            'benzimidazole': 'C12=CC=CC=C1NC=N2',
            'biphenyltetrazole': 'C1(C2=C(C3=CC=CC=C3)C=CC=C2)=NN=NN1',
            '3,3-hydroxy-2-oxindole': 'OC(C1=CC=CC=C1N2)C2=O',
            '5,7,5-lactone': 'C=C1C2CCCC3C(CC3)C2OC1=O',
            '6,6-spiroacetal': 'C1CCCC2(CCCCO2)O1',
            'dihydropyrimidone': 'O=C1NCC=CN1',
            'indolizine': 'N12C=CC=C1C=CC=C2',
            'biphenyl': 'C1(C2=CC=CC=C2)=CC=CC=C1',
            'triazaspirodecanone': 'O=C(NC1)C2(CCNCC2)N1C3=CC=CC=C3',
            'N-acylhydrazone': '[H]C(/N=N/CC)=O',
            'pyrrolinone': 'O=C1C=CNC1',
            'hydroxyamate': 'ONC(CCC(C)=O)=O',
            'trans-lactam': 'O=C1NC2CCCC2C1',
            'trans-lactone': 'O=C1OC2CCCC2C1',
            'hexahydroisoindole': 'C12CNCC1CCCC2',
            'benzimidazolone': 'O=C1N(C2CCNCC2)C3=CC=CC=C3N1',
            'indoline': 'C12=C(NCC2)C=CC=C1',
            '2-arylbenzothiazole': 'C12=CC=CC=C1N=C(C3=CC=CC=C3)S2',
            'imidazolequinoxaline': 'C1(NC2)=CC=CC=C1N3C2=CN=C3',
            'spiroindanylpiperidine': 'C12=CC=CC=C1C3(CCNCC3)CC2',
            'aminopyridazine': 'NC1=NN=CC=C1',
            '1,4-pyrazolodiazepin-8-one': 'O=C1NCCNC2=CNN=C21',
            'rhodanine': 'S=C(N1)SCC1=O',
            'pyranopyridone': 'O=C1C2=C(OCC=C2)C=CN1',
            'pyranoquinolone': 'O=C1C=CC2=CC=CC=C2N1'
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'indole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#7H]:2',
            'quinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6]:[#6]:[#6]:2',
            'isoquinoline': '[#6]12:[#6](:[#6]:[#7]:[#6]:[#6]:1):[#6]:[#6]:[#6]:[#6]:2',
            'purine': '[#6]12:[#7]:[#6]:[#7]:[#6]:[#6]:1:[#7H]:[#6]:[#7]:2',
            'quinoxaline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6]:[#6]:[#7]:2',
            'quinazolinone': '[#8]=[#6]1:[#7H]:[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#7]:1',
            'tetrahydroisoquinoline': '[#6]12:[#6](-[#6]-[#7]-[#6]-[#6]-1):[#6]:[#6]:[#6]:[#6]:2',
            'tetrahydraquinoline': '[#6]12:[#6](-[#7]-[#6]-[#6]-[#6]-1):[#6]:[#6]:[#6]:[#6]:2',
            'benzoxazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#8]:[#6]:[#7]:2',
            'benzofuran': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#8]:2',
            '3,3-dimethylbenzopyran': '[#6]-[#6]1(-[#6])-[#6]=[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#8]-1',
            'chromone': '[#8]=[#6]1:[#6]:[#6]:[#8]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
            'coumarin': '[#8]=[#6]1:[#8]:[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
            'carbohydrate': '[#8]-[#6]-[#6]1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]',
            'steroid': '[#6]12-[#6]-[#6]-[#6]-[#6]-[#6]-1-[#6]1-[#6](-[#6]3-[#6]-[#6]-[#6]-[#6]-3-[#6]-[#6]-1)-[#6]-[#6]-2',
            'prostanoic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6@@H]1-[#6@H](-[#6]-[#6]-[#6]-1)-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#8])=[#8]',
            'benzodiazepine': '[#8]=[#6]1-[#6]-[#7]=[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2)-[#7]-1',
            'arylpiperidine': '[#6]1(-[#6]2-[#6]-[#6]-[#7]-[#6]-[#6]-2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'arylpiperizine': '[#6]1(-[#7]2-[#6]-[#6]-[#7]-[#6]-[#6]-2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'benzylpiperidine': '[#7]1(-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            'benzothiophene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#16]:2',
            'dihydropyridine': '[#6]1-[#6]-[#6]=[#6]-[#6]=[#7]-1',
            'benzimidazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7H]:[#6]:[#7]:2',
            'biphenyltetrazole': '[#6]1(-[#6]2:[#6](-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2):[#7]:[#7]:[#7]:[#7H]:1',
            '3,3-hydroxy-2-oxindole': '[#8]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#7]-[#6]-1=[#8]',
            '5,7,5-lactone': '[#6]=[#6]1-[#6]2-[#6]-[#6]-[#6]-[#6]3-[#6](-[#6]-[#6]-3)-[#6]-2-[#8]-[#6]-1=[#8]',
            '6,6-spiroacetal': '[#6]1-[#6]-[#6]-[#6]-[#6]2(-[#6]-[#6]-[#6]-[#6]-[#8]-2)-[#8]-1',
            'dihydropyrimidone': '[#8]=[#6]1-[#7]-[#6]-[#6]=[#6]-[#7]-1',
            'indolizine': '[#7]12:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
            'biphenyl': '[#6]1(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'triazaspirodecanone': '[#8]=[#6]1-[#7]-[#6]-[#7](-[#6]-12-[#6]-[#6]-[#7]-[#6]-[#6]-2)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'N-acylhydrazone': '[#6H](/[#7]=[#7]/[#6]-[#6])=[#8]',
            'pyrrolinone': '[#8]=[#6]1-[#6]=[#6]-[#7]-[#6]-1',
            'hydroxyamate': '[#8]-[#7]-[#6](-[#6]-[#6]-[#6](-[#6])=[#8])=[#8]',
            'trans-lactam': '[#8]=[#6]1-[#7]-[#6]2-[#6]-[#6]-[#6]-[#6]-2-[#6]-1',
            'trans-lactone': '[#8]=[#6]1-[#8]-[#6]2-[#6]-[#6]-[#6]-[#6]-2-[#6]-1',
            'hexahydroisoindole': '[#6]12-[#6]-[#7]-[#6]-[#6]-1-[#6]-[#6]-[#6]-[#6]-2',
            'benzimidazolone': '[#8]=[#6]1:[#7](-[#6]2-[#6]-[#6]-[#7]-[#6]-[#6]-2):[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#7H]:1',
            'indoline': '[#6]12:[#6](-[#7]-[#6]-[#6]-1):[#6]:[#6]:[#6]:[#6]:2',
            '2-arylbenzothiazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1):[#16]:2',
            'imidazolequinoxaline': '[#6]12-[#7]-[#6]-[#6]3:[#7](-[#6]:1:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#7]:[#6]:3',
            'spiroindanylpiperidine': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]1(-[#6]-[#6]-[#7]-[#6]-[#6]-1)-[#6]-[#6]-2',
            'aminopyridazine': '[#7]-[#6]1:[#7]:[#7]:[#6]:[#6]:[#6]:1',
            '1,4-pyrazolodiazepin-8-one': '[#8]=[#6]1-[#7]-[#6]-[#6]-[#7]-[#6]2:[#6]:[#7H]:[#7]:[#6]:2-1',
            'rhodanine': '[#16]=[#6]1-[#7]-[#6](-[#6]-[#16]-1)=[#8]',
            'pyranopyridone': '[#8]=[#6]1:[#6]2:[#6](-[#8]-[#6]-[#6]=[#6]-2):[#6]:[#6]:[#7H]:1',
            'pyranoquinolone': '[#8]=[#6]1:[#6]:[#6]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#7H]:1',
        }

        return smarts
