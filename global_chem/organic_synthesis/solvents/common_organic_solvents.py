#!/usr/bin/env python3
#
# GlobalChem - Common Organic Solvents
#
# -----------------------------------

class CommonOrganicSolvents(object):

    def __init__(self):

        self.name = 'common_organic_solvents'

    @staticmethod
    def get_smiles():

        smiles = {
            'acetic acid': 'CC(=O)O',
            'acetone' : 'CC(=O)C',
            'acetonitrile': 'CC#N',
            'benzene': 'C1=CC=CC=C1',
            'tert-butyl alcohol': 'CC(C)(C)O',
            'tert-butyl methyl ether': 'CC(C)(C)OC',
            'butylated hydroxytoluene': 'CC1=CC(=C(C(=C1)C(C)(C)C)O)C(C)(C)C',
            'chloroform': 'C(Cl)(Cl)Cl',
            '18-crown-6': 'C1COCCOCCOCCOCCOCCO1',
            'cyclohexane': 'C1CCCCC1',
            '1,2-dichloroethane': 'C(CCl)Cl',
            'dichloromethane': 'C(Cl)Cl',
            'diethyl ether': 'CCOCC',
            'diglyme': 'COCCOCCOC',
            '1,2-dimethoxyethane': 'COCCOC',
            'dimethylacetamide': 'CC(=O)N(C)C',
            'dimethylformamide': 'CN(C)C=O',
            'dimethyl sulfoxide': 'CS(=O)C',
            'dioxane': 'C1COCCO1',
            'ethanol': 'CCO',
            'ethyl acetate': 'CCOC(=O)C',
            'ethyl methyl ketone': 'CCC(=O)C',
            'ethylene': 'C=C',
            'ethylene glycol': 'C(CO)O',
            'grease': 'C(C(F)(F)F)OCC(F)(F)F',
            'n-hexane': 'CCCCCC',
            'hexamethylbenzene': 'CC1=C(C(=C(C(=C1C)C)C)C)C',
            'hexamethylphosphoramide': 'CN(C)P(=O)(N(C)C)N(C)C',
            'hexamethyldisiloxane': 'O([Si](C)(C)C)[Si](C)(C)C',
            'methanol': 'CO',
            'nitromethane': 'C[N+](=O)[O-]',
            'n-pentane': 'CCCCC',
            'propylene': 'CC=C',
            '2-propanol': 'CC(C)O',
            'pyridine': 'C1=CC=NC=C1',
            'pyrrole': 'C1=CNC=C1',
            'pyrrolidine': 'C1CCNC1',
            'silicon grease': 'C[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C',
            'tetrahydrofuran': 'C1CCOC1',
            'toluene': 'CC1=CC=CC=C1',
            'triethylamine': 'CCN(CC)CC',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'acetic acid': '[#6]-[#6](=[#8])-[#8]',
            'acetone': '[#6]-[#6](=[#8])-[#6]',
            'acetonitrile': '[#6]-[#6]#[#7]',
            'benzene': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'tert-butyl alcohol': '[#6]-[#6](-[#6])(-[#6])-[#8]',
            'tert-butyl methyl ether': '[#6]-[#6](-[#6])(-[#6])-[#8]-[#6]',
            'butylated hydroxytoluene': '[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#6](-[#6])(-[#6])-[#6])-[#8])-[#6](-[#6])(-[#6])-[#6]',
            'chloroform': '[#6](-[#17])(-[#17])-[#17]',
            '18-crown-6': '[#6]1-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-1',
            'cyclohexane': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            '1,2-dichloroethane': '[#6](-[#6]-[#17])-[#17]',
            'dichloromethane': '[#6](-[#17])-[#17]',
            'diethyl ether': '[#6]-[#6]-[#8]-[#6]-[#6]',
            'diglyme': '[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]',
            '1,2-dimethoxyethane': '[#6]-[#8]-[#6]-[#6]-[#8]-[#6]',
            'dimethylacetamide': '[#6]-[#6](=[#8])-[#7](-[#6])-[#6]',
            'dimethylformamide': '[#6]-[#7](-[#6])-[#6]=[#8]',
            'dimethyl sulfoxide': '[#6]-[#16](=[#8])-[#6]',
            'dioxane': '[#6]1-[#6]-[#8]-[#6]-[#6]-[#8]-1',
            'ethanol': '[#6]-[#6]-[#8]',
            'ethyl acetate': '[#6]-[#6]-[#8]-[#6](=[#8])-[#6]',
            'ethyl methyl ketone': '[#6]-[#6]-[#6](=[#8])-[#6]',
            'ethylene': '[#6]=[#6]',
            'ethylene glycol': '[#6](-[#6]-[#8])-[#8]',
            'grease': '[#6](-[#6](-[#9])(-[#9])-[#9])-[#8]-[#6]-[#6](-[#9])(-[#9])-[#9]',
            'n-hexane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]',
            'hexamethylbenzene': '[#6]-[#6]1:[#6](:[#6](:[#6](:[#6](:[#6]:1-[#6])-[#6])-[#6])-[#6])-[#6]',
            'hexamethylphosphoramide': '[#6]-[#7](-[#6])-[#15](=[#8])(-[#7](-[#6])-[#6])-[#7](-[#6])-[#6]',
            'hexamethyldisiloxane': '[#8](-[Si](-[#6])(-[#6])-[#6])-[Si](-[#6])(-[#6])-[#6]',
            'methanol': '[#6]-[#8]',
            'nitromethane': '[#6]-[#7+](=[#8])-[#8-]',
            'n-pentane': '[#6]-[#6]-[#6]-[#6]-[#6]',
            'propylene': '[#6]-[#6]=[#6]',
            '2-propanol': '[#6]-[#6](-[#6])-[#8]',
            'pyridine': '[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1',
            'pyrrole': '[#6]1:[#6]:[#7H]:[#6]:[#6]:1',
            'pyrrolidine': '[#6]1-[#6]-[#6]-[#7]-[#6]-1',
            'silicon grease': '[#6]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#6]',
            'tetrahydrofuran': '[#6]1-[#6]-[#6]-[#8]-[#6]-1',
            'toluene': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'triethylamine': '[#6]-[#6]-[#7](-[#6]-[#6])-[#6]-[#6]'
        }

        return smarts
