#!/usr/bin/env python3
#
# GlobalChem - Privileged Kinase Inhibitor Scaffolds
#
# --------------------------------------------------

class PrivilegedKinaseInhibitors(object):

    def __init__(self):

        pass

    @staticmethod
    def get_smiles():

        smiles = {
            'indole': 'C12=CC=CC=C1C=CN2',
            'quinoline': 'C12=CC=CC=C1C=CC=N2',
            'phenylpiperazine': 'C1(N2CCNCC2)=CC=CC=C1',
            'biphenyl': 'C1(C2=CC=CC=C2)=CC=CC=C1',
            'benzimidazole': 'C12=CC=CC=C1NC=N2',
            'quinazoline': 'C12=CC=CC=C1C=NC=N2',
            'purine': 'C12=NC=NC=C1NC=N2',
            'indoline': 'C12=CC=CC=C1CCN2',
            'isoquinoline': 'C12=CC=CC=C1C=NC=C2',
            'benzylpiperidine': 'N1(CC2=CC=CC=C2)CCCCC1',
            'aminopyridazine': 'NC1=CC=CN=N1',
            '4-phenylpiperidine': 'C1(C2CCNCC2)=CC=CC=C1',
            'chromone': 'O=C1C=COC2=CC=CC=C21',
            '4-hydroxyquinazoline': 'O=C1NC=NC2=CC=CC=C21',
            'benzothiophene': 'C12=CC=CC=C1SC=C2',
            'benzofuran': 'C12=CC=CC=C1OC=C2',
            'quinoxaline': 'C12=CC=CC=C1N=CC=N2',
            'benzo[d]oxazole': 'C12=CC=CC=C1OC=N2',
            '1,2,3,4-tetrahydroisoquinoline': 'C12=CC=CC=C1CCNC2',
            'thiazolidine-2,4-dione': 'O=C(N1)SCC1=O',
            '1,2,3,4-tetrahydroquinoline': 'C12=CC=CC=C1CCCN2',
            '2H-chromen-2-one': 'O=C1OC2=CC=CC=C2C=C1',
            '1-(piperidin-4-yl)-1,3-dihydro-2H-benzo[d]imidazol-2-one': 'O=C1N(C2CCNCC2)C3=CC=CC=C3N1',
            '5H-dibenzo[b,e][1,4]diazepine': 'C1(C=CC=C2)=C2NC(C=CC=C3)=C3C=N1',
            '3,4-dihydropyrimidin-2(1H)-one': 'O=C1NC=CCN1',
            '3,4-dihydropyrimidine-2(1H)-thione': 'S=C1NC=CCN1',
            '6-(hydroxymethyl)tetrahydro-2H-pyran-2,3,4,5-tetraol': 'OCC1OC(O)C(O)C(O)C1O',
            '1-phenyl-1,3,8-triazaspiro[4.5]decan-4-one': 'O=C(NC1)C2(CCNCC2)N1C3=CC=CC=C3',
            '1,4-dihydropyridine': 'C1=CNC=CC1',
            '2-(tetrazol-5-yl)biphenyl': 'C1(C2=CC=CC=C2C3=NN=NN3)=CC=CC=C1'
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'indole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#7H]:2',
            'quinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#6]:[#7]:2',
            'phenylpiperazine': '[#6]1(-[#7]2-[#6]-[#6]-[#7]-[#6]-[#6]-2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'biphenyl': '[#6]1(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'benzimidazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7H]:[#6]:[#7]:2',
            'quinazoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#6]:[#7]:2',
            'purine': '[#6]12:[#7]:[#6]:[#7]:[#6]:[#6]:1:[#7H]:[#6]:[#7]:2',
            'indoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#7]-2',
            'isoquinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#6]:[#6]:2',
            'benzylpiperidine': '[#7]1(-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            'aminopyridazine': '[#7]-[#6]1:[#6]:[#6]:[#6]:[#7]:[#7]:1',
            '4-phenylpiperidine': '[#6]1(-[#6]2-[#6]-[#6]-[#7]-[#6]-[#6]-2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'chromone': '[#8]=[#6]1:[#6]:[#6]:[#8]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:1:2',
            '4-hydroxyquinazoline': '[#8]=[#6]1:[#7H]:[#6]:[#7]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:1:2',
            'benzothiophene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#16]:[#6]:[#6]:2',
            'benzofuran': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#8]:[#6]:[#6]:2',
            'quinoxaline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6]:[#6]:[#7]:2',
            'benzo[d]oxazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#8]:[#6]:[#7]:2',
            '1,2,3,4-tetrahydroisoquinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#7]-[#6]-2',
            'thiazolidine-2,4-dione': '[#8]=[#6]1-[#7]-[#6](-[#6]-[#16]-1)=[#8]',
            '1,2,3,4-tetrahydroquinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]-[#7]-2',
            '2H-chromen-2-one': '[#8]=[#6]1:[#8]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#6]:[#6]:1',
            '1-(piperidin-4-yl)-1,3-dihydro-2H-benzo[d]imidazol-2-one': '[#8]=[#6]1:[#7](-[#6]2-[#6]-[#6]-[#7]-[#6]-[#6]-2):[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#7H]:1',
            '5H-dibenzo[b,e][1,4]diazepine': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]=[#7]-2',
            '3,4-dihydropyrimidin-2(1H)-one': '[#8]=[#6]1-[#7]-[#6]=[#6]-[#6]-[#7]-1',
            '3,4-dihydropyrimidine-2(1H)-thione': '[#16]=[#6]1-[#7]-[#6]=[#6]-[#6]-[#7]-1',
            '6-(hydroxymethyl)tetrahydro-2H-pyran-2,3,4,5-tetraol': '[#8]-[#6]-[#6]1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]',
            '1-phenyl-1,3,8-triazaspiro[4.5]decan-4-one': '[#8]=[#6]1-[#7]-[#6]-[#7](-[#6]-12-[#6]-[#6]-[#7]-[#6]-[#6]-2)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '1,4-dihydropyridine': '[#6]1=[#6]-[#7]-[#6]=[#6]-[#6]-1',
            '2-(tetrazol-5-yl)biphenyl': '[#6]1(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]2:[#7]:[#7]:[#7]:[#7H]:2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
        }

        return smarts
