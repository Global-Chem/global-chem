#!/usr/bin/env python3
#
# GlobalChem - Cimetidine & Acyclovir
#
# -----------------------------------

class CimetidineAndAcyclovir(object):

    def __init__(self):

        self.name = 'cimetidine_and_acyclovir'

    @staticmethod
    def get_smiles():

        smiles = {
            'microcrystalline cellulose': 'COC1OC(CO)C(OC2OC(CO)C(OC)C(O)C2O)C(O)C1O',
            'corn starch': 'COC4C(CO)OC(OCC2OC(OC1C(CO)OC(C)C(O)C1O)C(O)C(O)C2OC3OC(CO)C(OC)C(O)C3O)C(O)C4O',
            'dibasic calcium phosphate': 'O.O.OP(=O)([O-])[O-].[Ca+2]',
            'lactose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
            'pregelatinized starch': 'COC1C(O)C(O)C(OCC2OC(OC3C(O)C(O)C(C)OC3CO)C(O)C(O)C2OC2OC(CO)C(OC)C(O)C2O)OC1CO',
            'hydroxypropyl methylcellulose': 'CC(COCC1C(C(C(C(O1)OC2C(OC(C(C2OCC(C)O)OCC(C)O)OCC(C)O)COCC(C)O)OCC(C)O)OCC(C)O)OCC(C)O)O.COCC1C(C(C(C(O1)OC2C(OC(C(C2OC)OC)OC)COC)OC)OC)OC',
            'sodium starch glycolate': 'NN1CC2N(C(C)=O)C(C1)CC2',
            'sodium lauryl sulfate': 'CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]',
            'povidone': 'C1CC(=O)N(C1)C(CP)P',
            'croscarmellose_sodium': 'CC(=O)O.C(C(C(C(C(C=O)O)O)O)O)O.[Na]',
            'colloidal_silicon_dioxide': 'O=[Si]=O',
            'crospovidone': 'C1CC(=O)NC1',
            'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O',
            'magnesium stearate': 'CCCCCCCCCCCCCCCCCC(=O)[O-].CCCCCCCCCCCCCCCCCC(=O)[O-].[Mg+2]'
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'microcrystalline cellulose': '[#6]-[#8]-[#6]1-[#8]-[#6](-[#6]-[#8])-[#6](-[#8]-[#6]2-[#8]-[#6](-[#6]-[#8])-[#6](-[#8]-[#6])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])-[#6]-1-[#8]',
            'corn starch': '[#6]-[#8]-[#6]1-[#6](-[#6]-[#8])-[#8]-[#6](-[#8]-[#6]-[#6]2-[#8]-[#6](-[#8]-[#6]3-[#6](-[#6]-[#8])-[#8]-[#6](-[#6])-[#6](-[#8])-[#6]-3-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8]-[#6]2-[#8]-[#6](-[#6]-[#8])-[#6](-[#8]-[#6])-[#6](-[#8])-[#6]-2-[#8])-[#6](-[#8])-[#6]-1-[#8]',
            'dibasic calcium phosphate': '[#8].[#8].[#8]-[#15](=[#8])(-[#8-])-[#8-].[Ca+2]',
            'lactose': '[#6](-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1-[#6](-[#8]-[#6](-[#6](-[#6]-1-[#8])-[#8])-[#8])-[#6]-[#8])-[#8])-[#8])-[#8])-[#8]',
            'pregelatinized starch': '[#6]-[#8]-[#6]1-[#6](-[#8])-[#6](-[#8])-[#6](-[#8]-[#6]-[#6]2-[#8]-[#6](-[#8]-[#6]3-[#6](-[#8])-[#6](-[#8])-[#6](-[#6])-[#8]-[#6]-3-[#6]-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-2-[#8]-[#6]2-[#8]-[#6](-[#6]-[#8])-[#6](-[#8]-[#6])-[#6](-[#8])-[#6]-2-[#8])-[#8]-[#6]-1-[#6]-[#8]',
            'hydroxypropyl methylcellulose': '[#6]-[#6](-[#6]-[#8]-[#6]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1-[#6](-[#8]-[#6](-[#6](-[#6]-1-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#6]-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8].[#6]-[#8]-[#6]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1-[#6](-[#8]-[#6](-[#6](-[#6]-1-[#8]-[#6])-[#8]-[#6])-[#8]-[#6])-[#6]-[#8]-[#6])-[#8]-[#6])-[#8]-[#6])-[#8]-[#6]',
            'sodium starch glycolate': '[#7]-[#7]1-[#6]-[#6]2-[#7](-[#6](-[#6])=[#8])-[#6](-[#6]-1)-[#6]-[#6]-2',
            'sodium lauryl sulfate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#8]-[#16](=[#8])(=[#8])-[#8-].[Na+]',
            'povidone': '[#6]1-[#6]-[#6](=[#8])-[#7](-[#6]-1)-[#6](-[#6]-[#15])-[#15]',
            'croscarmellose sodium': '[#6]-[#6](=[#8])-[#8].[#6](-[#6](-[#6](-[#6](-[#6](-[#6]=[#8])-[#8])-[#8])-[#8])-[#8])-[#8].[Na]',
            'colloidal silicon dioxide': '[#8]=[Si]=[#8]',
            'crospovidone': '[#6]1-[#6]-[#6](=[#8])-[#7]-[#6]-1',
            'stearic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]',
            'magnesium stearate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8-].[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8-].[Mg+2]',
        }

        return smarts