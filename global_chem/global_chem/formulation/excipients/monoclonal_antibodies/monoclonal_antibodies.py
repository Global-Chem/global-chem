#!/usr/bin/env python3
#
# GlobalChem - Monoclonal Antibodies
#
# -----------------------------------

class MonoclonalAntibodies(object):

    def __init__(self):

        self.name = 'monoclonal_antibodies'

    @staticmethod
    def get_smiles():

        smiles = {
            'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
            'polysorbate 80': r'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
            'sucrose': 'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O',
            'glycine': 'C(C(=O)O)N',
            'sodium phosphate, dibasic, anhydrous': 'OP(=O)([O-])[O-].[Na+].[Na+]',
            'potassium phosphate monobasic': 'OP(=O)(O)[O-].[K+]',
            'sodium phosphate, monobasic, monohydrate': 'OP(=O)(O)[O-].[Na+]',
            'sodium phosphate, dibasic, dihydrate': 'O.O.OP(=O)([O-])[O-].[Na+].[Na+]',
            'dibasic sodium phosphate heptahydrate': 'O.O.O.O.O.O.O.OP(=O)([O-])[O-].[Na+].[Na+] ',
            'monobasic sodium phosphate monohydrate': 'O.OP(=O)(O)[O-].[Na+]',
            'histidine': 'C1=C(NC=N1)CC(C(=O)O)N',
            'sodium citrate dihydrate': 'C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O.O.O.[Na+].[Na+].[Na+]',
            'histidine hydrochloride monohydrate': 'C1=C(NC=N1)CC(C(=O)O)N.O.Cl',
            'trehalose': 'C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O',
            'edetate disodium dihydrate': 'C(CN(CC(=O)O)CC(=O)[O-])N(CC(=O)O)CC(=O)[O-].O.O.[Na+].[Na+]',
            'sorbitol': 'C(C(C(C(C(CO)O)O)O)O)O',
            'polysorbate 20': 'CCCCCCCCCCCC(=O)OCCOCC(C1C(C(CO1)OCCO)OCCO)OCCO',
            'succinic acid': 'C(CC(=O)O)C(=O)O',
            'ethylene glycol': 'C(CO)O',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'mannitol': '[#6](-[#6](-[#6](-[#6](-[#6](-[#6]-[#8])-[#8])-[#8])-[#8])-[#8])-[#8]',
            'polysorbate 80': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]/[#6]=[#6]/[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6](-[#6]1-[#6](-[#6](-[#6]-[#8]-1)-[#8])-[#8])-[#8]',
            'sucrose': '[#6](-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1(-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#6]-[#8])-[#8])-[#8])-[#8])-[#8]',
            'glycine': '[#6](-[#6](=[#8])-[#8])-[#7]',
            'sodium phosphate, dibasic, anhydrous': '[#8]-[#15](=[#8])(-[#8-])-[#8-].[Na+].[Na+]',
            'potassium phosphate monobasic': '[#8]-[#15](=[#8])(-[#8])-[#8-].[K+]',
            'sodium phosphate, monobasic, monohydrate': '[#8]-[#15](=[#8])(-[#8])-[#8-].[Na+]',
            'sodium phosphate, dibasic, dihydrate': '[#8].[#8].[#8]-[#15](=[#8])(-[#8-])-[#8-].[Na+].[Na+]',
            'dibasic sodium phosphate heptahydrate': '[#8].[#8].[#8].[#8].[#8].[#8].[#8].[#8]-[#15](=[#8])(-[#8-])-[#8-].[Na+].[Na+]',
            'monobasic sodium phosphate monohydrate': '[#8].[#8]-[#15](=[#8])(-[#8])-[#8-].[Na+]',
            'histidine': '[#6]1:[#6](:[#7H]:[#6]:[#7]:1)-[#6]-[#6](-[#6](=[#8])-[#8])-[#7]',
            'sodium citrate dihydrate': '[#6](-[#6](=[#8])-[#8-])-[#6](-[#6]-[#6](=[#8])-[#8-])(-[#6](=[#8])-[#8-])-[#8].[#8].[#8].[Na+].[Na+].[Na+]',
            'histidine hydrochloride monohydrate': '[#6]1:[#6](:[#7H]:[#6]:[#7]:1)-[#6]-[#6](-[#6](=[#8])-[#8])-[#7].[#8].[#17]',
            'trehalose': '[#6](-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8])-[#8])-[#8])-[#8])-[#8]',
            'edetate disodium dihydrate': '[#6](-[#6]-[#7](-[#6]-[#6](=[#8])-[#8])-[#6]-[#6](=[#8])-[#8-])-[#7](-[#6]-[#6](=[#8])-[#8])-[#6]-[#6](=[#8])-[#8-].[#8].[#8].[Na+].[Na+]',
            'sorbitol': '[#6](-[#6](-[#6](-[#6](-[#6](-[#6]-[#8])-[#8])-[#8])-[#8])-[#8])-[#8]',
            'polysorbate 20': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6]-[#8]-[#6]-[#6](-[#6]1-[#6](-[#6](-[#6]-[#8]-1)-[#8]-[#6]-[#6]-[#8])-[#8]-[#6]-[#6]-[#8])-[#8]-[#6]-[#6]-[#8]',
            'succinic acid': '[#6](-[#6]-[#6](=[#8])-[#8])-[#6](=[#8])-[#8]',
            'ethylene glycol': '[#6](-[#6]-[#8])-[#8]',
        }

        return smarts