#!/usr/bin/env python3
#
# GlobalChem - Braf Inhibitor by Pocket
#
# -------------------------------------

class BRAFInhibitors(object):

    def __init__(self):

        pass

    @staticmethod
    def get_smiles():

        smiles = {
            'imidazole': 'C1=CN=CN1',
            'pyrazole-1-ethanol': 'OCCN1N=CC=C1',
            'acetonitrile': 'C#N',
            '2-(tert-butyl)thiazole': 'CC(C)(C)C1=NC=CS1',
            'pyridine': 'C1=CC=NC=C1',
            '1-isopropyl-pyrazole': 'CC(C)N1N=CC=C1',
            'isoindoline' :'C12=CC=CC=C1CNC2',
            'pyrrolopyridine': 'C12=CC=CN=C1C=CN2',
            '1H-pyrrolo[2,3-b]pyridine': 'C12=NC=CC=C1C=CN2',
            'pyrimidine': 'C1=CN=CN=C1',
            '3,4-dihydroquinazoline': 'C12=CC=CC=C1CNC=N2',
            '1,3,8-Triazanaphthalene': 'C12=NC=CC=C1C=NC=N2',
            'benzothiazole': 'C1(SC=N2)=C2C=CC=C1',
            '3,4-dihydropyrido[2,3-b]pyrazine':'C12=NC=CC=C1N=CCN2',
            'morpholine': 'C1COCCN1',
            '2-aminopyrimidine': 'NC1=NC=CC=N1',
            'benzoimidazole': 'C12=CC=CC=C1N=CN2',
            'pyridinyl imidazole': 'C1(C2=NC=CN2)=CC=CC=N1',
            'quinazoline': 'C12=NC=NC=C1C=CC=C2',
            'benzene': 'C1=CC=CC=C1',
            '1-methyleneindan': 'C=C1CCC2=C1C=CC=C2',
            'toluene': 'CC1=CC=CC=C1',
            'fluorobenzene': 'FC1=CC=CC=C1',
            'methyl(phenyl)sulfane': 'CSC1=CC=CC=C1',
            '2-methylpyridine': 'CC1=CC=CC=N1',
            '1a,6b-dihydro-1H-cyclopropa[b]benzofuran': 'C12=CC=CC=C1OC3C2C3',
            '1,2,3,4-tetrahydronaphthalene': 'C12=CC=CC=C1CCCC2',
            '1,3-difluorobenzene': 'FC1=CC(F)=CC=C1',
            '1-chloro-4-fluorobenzene': 'FC1=CC=C(Cl)C=C1',
            'sulfonylpropane': 'O=S(CCC)=O',
            'methanesulfonamide': 'O=S(C)(N)=O',
            '1,3-difluoro-2-sulfonylbenzene' : 'O=S(C1=C(F)C=CC=C1F)=O',
            'N-ethyl-N-methylsulfonic amide': 'O=S(N(CC)C)=O',
            'propanesulfonamide': 'O=S(CCC)(N)=O',
            '1-hydrosulfonylpyrrolidine': 'O=S(N1CCCC1)=O',
            'prop-2-yn-1-ylbenzene': 'C#CCC1=CC=CC=C1',
            'neohexane': 'CCC(C)(C)C',
            '(trifluoromethyl)benzene': 'FC(F)(F)C1=CC=CC=C1',
            '3-chloro-4-(trifluoromethyl)pyridine': 'ClC1=CN=CC=C1C(F)(F)F',
            '1-phenyl-1H-pyrazole': 'C1(N2N=CC=C2)=CC=CC=C1',
            '4-(trifluoromethyl)pyridine': 'FC(F)(F)C1=CC=NC=C1',
            '2-phenylpropan-2-amine': 'CC(C)(N)C1=CC=CC=C1',
            '6-methyl-1H-benzoimidazole': 'CC1=CC=C(N=CN2)C2=C1',
            '5-(1,1,1-trifluoro-2-methylpropan-2-yl)isoxazole': 'CC(C(F)(F)F)(C)C1=CC=NO1',
            'N,N-dimethyl-2-phenoxyethan-1-amine': 'CN(C)CCOC1=CC=CC=C1',
            'cyclopropanecarboxamide': 'O=C(N)C1CC1',
            'chlorobenzene': 'ClC1=CC=CC=C1',
            'methyl tert-butylcarbamate': 'COC(NC(C)(C)C)=O',
            'tetrahydro-2H-pyran': 'C1CCCCO1',
            '2-cyclopropylpyrimidine': 'C1(C2CC2)=NC=CC=N1',
            '1-ethylpiperidine': 'CCN1CCCCC1'
        }

        return smiles

    @staticmethod
    def get_smarts():
        
        smarts = {
            'imidazole': '[#6]1:[#6]:[#7]:[#6]:[#7H]:1',
            'pyrazole-1-ethanol': '[#8]-[#6]-[#6]-[#7]1:[#7]:[#6]:[#6]:[#6]:1',
            'acetonitrile': '[#6]#[#7]',
            '2-(tert-butyl)thiazole': '[#6]-[#6](-[#6])(-[#6])-[#6]1:[#7]:[#6]:[#6]:[#16]:1',
            'pyridine': '[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1',
            '1-isopropyl-pyrazole': '[#6]-[#6](-[#6])-[#7]1:[#7]:[#6]:[#6]:[#6]:1',
            'isoindoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#7]-[#6]-2',
            'pyrrolopyridine': '[#6]12:[#6]:[#6]:[#6]:[#7]:[#6]:1:[#6]:[#6]:[#7H]:2',
            '1H-pyrrolo[2,3-b]pyridine': '[#6]12:[#7]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#7H]:2',
            'pyrimidine': '[#6]1:[#6]:[#7]:[#6]:[#7]:[#6]:1',
            '3,4-dihydroquinazoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#7]-[#6]=[#7]-2',
            '1,3,8-Triazanaphthalene': '[#6]12:[#7]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#6]:[#7]:2',
            'benzothiazole': '[#6]12:[#16]:[#6]:[#7]:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
            '3,4-dihydropyrido[2,3-b]pyrazine': '[#6]12:[#7]:[#6]:[#6]:[#6]:[#6]:1-[#7]=[#6]-[#6]-[#7]-2',
            'morpholine': '[#6]1-[#6]-[#8]-[#6]-[#6]-[#7]-1',
            '2-aminopyrimidine': '[#7]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#7]:1',
            'benzoimidazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6]:[#7H]:2',
            'pyridinyl imidazole': '[#6]1(-[#6]2:[#7]:[#6]:[#6]:[#7H]:2):[#6]:[#6]:[#6]:[#6]:[#7]:1',
            'quinazoline': '[#6]12:[#7]:[#6]:[#7]:[#6]:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
            'benzene': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '1-methyleneindan': '[#6]=[#6]1-[#6]-[#6]-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2',
            'toluene': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'fluorobenzene': '[#9]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'methyl(phenyl)sulfane': '[#6]-[#16]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '2-methylpyridine': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1',
            '1a,6b-dihydro-1H-cyclopropa[b]benzofuran': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]-[#6]1-[#6]-2-[#6]-1',
            '1,2,3,4-tetrahydronaphthalene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]-[#6]-2',
            '1,3-difluorobenzene': '[#9]-[#6]1:[#6]:[#6](-[#9]):[#6]:[#6]:[#6]:1',
            '1-chloro-4-fluorobenzene': '[#9]-[#6]1:[#6]:[#6]:[#6](-[#17]):[#6]:[#6]:1',
            'sulfonylpropane': '[#8]=[#16](-[#6]-[#6]-[#6])=[#8]',
            'methanesulfonamide': '[#8]=[#16](-[#6])(-[#7])=[#8]',
            '1,3-difluoro-2-sulfonylbenzene': '[#8]=[#16](-[#6]1:[#6](-[#9]):[#6]:[#6]:[#6]:[#6]:1-[#9])=[#8]',
            'N-ethyl-N-methylsulfonic amide': '[#8]=[#16](-[#7](-[#6]-[#6])-[#6])=[#8]',
            'propanesulfonamide': '[#8]=[#16](-[#6]-[#6]-[#6])(-[#7])=[#8]',
            '1-hydrosulfonylpyrrolidine': '[#8]=[#16](-[#7]1-[#6]-[#6]-[#6]-[#6]-1)=[#8]',
            'prop-2-yn-1-ylbenzene': '[#6]#[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'neohexane': '[#6]-[#6]-[#6](-[#6])(-[#6])-[#6]',
            '(trifluoromethyl)benzene': '[#9]-[#6](-[#9])(-[#9])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '3-chloro-4-(trifluoromethyl)pyridine': '[#17]-[#6]1:[#6]:[#7]:[#6]:[#6]:[#6]:1-[#6](-[#9])(-[#9])-[#9]',
            '1-phenyl-1H-pyrazole': '[#6]1(-[#7]2:[#7]:[#6]:[#6]:[#6]:2):[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '4-(trifluoromethyl)pyridine': '[#9]-[#6](-[#9])(-[#9])-[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1',
            '2-phenylpropan-2-amine': '[#6]-[#6](-[#6])(-[#7])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '6-methyl-1H-benzoimidazole': '[#6]-[#6]1:[#6]:[#6]:[#6]2:[#7]:[#6]:[#7H]:[#6]:2:[#6]:1',
            '5-(1,1,1-trifluoro-2-methylpropan-2-yl)isoxazole': '[#6]-[#6](-[#6](-[#9])(-[#9])-[#9])(-[#6])-[#6]1:[#6]:[#6]:[#7]:[#8]:1',
            'N,N-dimethyl-2-phenoxyethan-1-amine': '[#6]-[#7](-[#6])-[#6]-[#6]-[#8]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'cyclopropanecarboxamide': '[#8]=[#6](-[#7])-[#6]1-[#6]-[#6]-1',
            'chlorobenzene': '[#17]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'methyl tert-butylcarbamate': '[#6]-[#8]-[#6](-[#7]-[#6](-[#6])(-[#6])-[#6])=[#8]',
            'tetrahydro-2H-pyran': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#8]-1',
            '2-cyclopropylpyrimidine': '[#6]1(-[#6]2-[#6]-[#6]-2):[#7]:[#6]:[#6]:[#6]:[#7]:1',
            '1-ethylpiperidine': '[#6]-[#6]-[#7]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        }
        
        return smarts
