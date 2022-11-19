#!/usr/bin/env python3
#
# GlobalChem - Alternative Jet Fuels
#
# ----------------------------------

class AlternativeJetFuels(object):

    def __init__(self):

        self.name = 'alternative_jet_fuels'

    @staticmethod
    def get_smiles():

        '''

        Missing Entries:

        '''

        smiles = {
            'heptane': 'CCCCCCC',
            '2-methyl heptane': 'CCCCCC(C)C',
            '3-methyl heptane': 'CCCCC(C)CC',
            'octane': 'CCCCCCCC',
            'ethyl cyclohexane': 'CCC1CCCCC1',
            '4-methyl octane': 'CCCCC(C)CCC',
            '2-methyl octane': 'CCCCCCC(C)C',
            '2,5-dimethyl heptane': 'CCC(C)CCC(C)C',
            '3-methyl octane': 'CCCCCC(C)CC',
            '2,2,4-trimethyl heptane': 'CCCC(C)CC(C)(C)C',
            'butyl cyclopentane': 'CCCCC1CCCC1',
            'nonane': 'CCCCCCCCC',
            '3-ethyl-2,5-dimethyl hexane': 'CCC(CC(C)C)C(C)C',
            '1,3-dimethyl benzene': 'CC1=CC(=CC=C1)C',
            'propyl cyclohexane': 'CCCC1CCCCC1',
            '2,5,6-trimethyl octane': 'CCC(C)C(C)CCC(C)C',
            '2-methyl nonane': 'CCCCCCCC(C)C',
            '3-methyl nonane': 'CCCCCCC(C)CC',
            '2,2,4,6,6-pentamethyl heptane': 'CC(CC(C)(C)C)CC(C)(C)C',
            '2,3-dimethyl nonane': 'CCCCCCC(C)C(C)C',
            'decane': 'CCCCCCCCCC',
            '1-methyl-2-propyl cyclohexane': 'CCCC1CCCCC1C',
            '5-ethyl-2,2,3-trimethyl heptane': 'CCC(CC)CC(C)C(C)(C)C',
            '1-ethyl-2-methyl benzene': 'CCC1=CC=CC=C1C',
            '3-methyl decane': 'CCCCCCCC(C)CC',
            '1,2,3-trimethyl benzene': 'CC1=C(C(=CC=C1)C)C',
            '3,7-dimethyl decane': 'CCCC(C)CCCC(C)CC',
            '4-methyl undecane': 'CCCCCCCC(C)CCC',
            'undecane': 'CCCCCCCCCCC',
            '1-ethenyl-2-methyl benzene': 'C=CC1=CC=CC=C1C',
            '1,2-diethyl benzene': 'CCC1=CC=CC=C1CC',
            '1-methyl-4-(2-propenyl) benzene': 'CC1=CC=C(C=C1)C=C(C)C',
            '1-ethyl-2,4-dimethyl benzene': 'CCC1=C(C=C(C=C1)C)C',
            '1-methyl-2-(4-methylpentyl)cyclopentane': 'CC1CCCC1CCCC(C)C',
            '3-methyl undecane': 'CCCCCCCCC(C)CC',
            '1,2,3,5-tetramethyl benzene': 'CC1=CC(=C(C(=C1)C)C)C',
            '2-methyl 1,1*-bicyclohexyl': 'CC1CCCCC1C2CCCCC2',
            '1,6-tridecadiene': 'CCCCCCC=CCCCC=C',
            'dodecane': 'CCCCCCCCCCCC',
            '1-ethenyl-4-ethyl benzene': 'CCC1=CC=C(C=C1)C=C',
            '1-methyl-2-(2-propenyl) benzene': 'CC=CC1=CC=CC=C1C',
            '2,6-dimethyl undecane': 'CCCCCC(C)CCCC(C)C',
            '1,2,3,4-tetrahydro naphthalene': 'C1CCC2=CC=CC=C2C1',
            '2-ethenyl-1,3,5-trimethyl benzene': 'C=CC1=C(C)C=C(C)C=C1C',
            'hexyl cyclohexane': 'CCCCCCC1CCCCC1',
            '3-methyl dodecane': 'CCCCCCCCCC(C)CC',
            '1,2,3,4-tetrahydro-2-methyl naphthalene': 'CC1CCC2=CC=CC=C2C1',
            'tridecane': 'CCCCCCCCCCCCC',
            '2,2,4,4,6,8,8-heptamethyl nonane': 'CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C',
            '1,2,3,4-tetrahydro-5-methyl naphthalene': 'CC1=C2CCCCC2=CC=C1',
            '3-methyl tridecane': 'CCCCCCCCCCC(C)CC',
            '2,6,10-trimethyl dodecane': 'CCC(C)CCCC(C)CCCC(C)C',
            '1,2,3,4-tetrahydro-2,7-diemthyl napthalene': 'CC1CCC2=C(C1)C=C(C=C2)C',
            'tetradecane': 'CCCCCCCCCCCCCC',
            '7-methyl pentadecane': 'CCCCCCCCC(C)CCCCCC',
            '3-methyl tetradecane': 'CCCCCCCCCCCC(C)CC',
            'pentadecane': 'CCCCCCCCCCCCCCC',
            'hexadecane': 'CCCCCCCCCCCCCCCC',
            'heptadecane': 'CCCCCCCCCCCCCCCCC',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'heptane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '2-methyl heptane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '3-methyl heptane': '[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            'octane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            'ethyl cyclohexane': '[#6]-[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            '4-methyl octane': '[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]',
            '2-methyl octane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '2,5-dimethyl heptane': '[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]',
            '3-methyl octane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '2,2,4-trimethyl heptane': '[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6](-[#6])(-[#6])-[#6]',
            'butyl cyclopentane': '[#6]-[#6]-[#6]-[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-1',
            'nonane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '3-ethyl-2,5-dimethyl hexane': '[#6]-[#6]-[#6](-[#6]-[#6](-[#6])-[#6])-[#6](-[#6])-[#6]',
            '1,3-dimethyl benzene': '[#6]-[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1)-[#6]',
            'propyl cyclohexane': '[#6]-[#6]-[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            '2,5,6-trimethyl octane': '[#6]-[#6]-[#6](-[#6])-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]',
            '2-methyl nonane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '3-methyl nonane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '2,2,4,6,6-pentamethyl heptane': '[#6]-[#6](-[#6]-[#6](-[#6])(-[#6])-[#6])-[#6]-[#6](-[#6])(-[#6])-[#6]',
            '2,3-dimethyl nonane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6](-[#6])-[#6]',
            'decane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '1-methyl-2-propyl cyclohexane': '[#6]-[#6]-[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1-[#6]',
            '5-ethyl-2,2,3-trimethyl heptane': '[#6]-[#6]-[#6](-[#6]-[#6])-[#6]-[#6](-[#6])-[#6](-[#6])(-[#6])-[#6]',
            '1-ethyl-2-methyl benzene': '[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]',
            '3-methyl decane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '1,2,3-trimethyl benzene': '[#6]-[#6]1:[#6](:[#6](:[#6]:[#6]:[#6]:1)-[#6])-[#6]',
            '3,7-dimethyl decane': '[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '4-methyl undecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]',
            'undecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '1-ethenyl-2-methyl benzene': '[#6]=[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]',
            '1,2-diethyl benzene': '[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]',
            '1-methyl-4-(2-propenyl) benzene': '[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]=[#6](-[#6])-[#6]',
            '1-ethyl-2,4-dimethyl benzene': '[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1)-[#6])-[#6]',
            '1-methyl-2-(4-methylpentyl)cyclopentane': '[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-1-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '3-methyl undecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '1,2,3,5-tetramethyl benzene': '[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#6])-[#6])-[#6]',
            '2-methyl 1,1*-bicyclohexyl': '[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            '1,6-tridecadiene': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]=[#6]-[#6]-[#6]-[#6]-[#6]=[#6]',
            'dodecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '1-ethenyl-4-ethyl benzene': '[#6]-[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]=[#6]',
            '1-methyl-2-(2-propenyl) benzene': '[#6]-[#6]=[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]',
            '2,6-dimethyl undecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '1,2,3,4-tetrahydro naphthalene': '[#6]1-[#6]-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]-1',
            '2-ethenyl-1,3,5-trimethyl benzene': '[#6]=[#6]-[#6]1:[#6](-[#6]):[#6]:[#6](-[#6]):[#6]:[#6]:1-[#6]',
            'hexyl cyclohexane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            '3-methyl dodecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '1,2,3,4-tetrahydro-2-methyl naphthalene': '[#6]-[#6]1-[#6]-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]-1',
            'tridecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '2,2,4,4,6,8,8-heptamethyl nonane': '[#6]-[#6](-[#6]-[#6](-[#6])(-[#6])-[#6])-[#6]-[#6](-[#6])(-[#6])-[#6]-[#6](-[#6])(-[#6])-[#6]',
            '1,2,3,4-tetrahydro-5-methyl naphthalene': '[#6]-[#6]1:[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]:2:[#6]:[#6]:[#6]:1',
            '3-methyl tridecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            '2,6,10-trimethyl dodecane': '[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]',
            '1,2,3,4-tetrahydro-2,7-diemthyl napthalene': '[#6]-[#6]1-[#6]-[#6]-[#6]2:[#6](-[#6]-1):[#6]:[#6](:[#6]:[#6]:2)-[#6]',
            'tetradecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '7-methyl pentadecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '3-methyl tetradecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
            'pentadecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            'hexadecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            'heptadecane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
        }

        return smarts
