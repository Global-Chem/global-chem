#!/usr/bin/env python3
#
# GlobalChem - Common Salts
#
# -----------------------------------

class Salt(object):

    def __init__(self):

        self.name = 'salt'

    @staticmethod
    def get_smiles():

        '''

        Missing:

            'concentré salin liquide': '',
            'eau et sel': '',
            'mineraalzout': '',
        '''

        smiles = {
            'siedesalz': '[Na+].[Cl-]',
            'black salt': '[Na+].[Cl-].[Na+].[O-]S([O-])(=O)=O.OS(=O)(=O)[O-].[O-]S(=O)(=O)[O-].[Fe+3]',
            'finegrained salt': '[Na+].[Cl-]',
            'himalayan salt': '[Na+].[Cl-].[K+].[Fe+3]',
            'iodised salt': '[Na+].[Cl-].[I-]',
            'kosher salt': '[Na+].[Cl-]',
            'pink salt': '[Na+].[Cl-].[K+].[Fe+3]',
            'potassium salt': '[K+].[Cl-].[I-]',
            'refined salt': '[Na+].[Cl-].[O-][Si](=O)[O-].[O-][Si](=O)[O-].[Al+3].C(=O)([O-])[O-].[Mg+2]',
            'sea salt': '[Na+].[K+].[Mg+2].[Ca+2].[Cl-]',
            'smoked salt': '[a].[Na+].[K+].[Mg+2].[Ca+2].[Cl-]',
            'sodium iodide': '[Na+].[I-]',
            'unrefined salt': '[Na+].[K+].[Mg+2].[Ca+2].[Cl-].[Zn+2].[Cu+].[Fe+3].[P]',
            'rock salt': '[Na+].[Cl-]',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'siedesalz': '[Na+].[#17-]',
            'black salt': '[Na+].[#17-].[Na+].[#8-]-[#16](-[#8-])(=[#8])=[#8].[#8]-[#16](=[#8])(=[#8])-[#8-].[#8-]-[#16](=[#8])(=[#8])-[#8-].[Fe+3]',
            'finegrained salt': '[Na+].[#17-]',
            'himalayan salt': '[Na+].[#17-].[K+].[Fe+3]',
            'iodised salt': '[Na+].[#17-].[#53-]',
            'kosher salt': '[Na+].[#17-]',
            'pink salt': '[Na+].[#17-].[K+].[Fe+3]',
            'potassium salt': '[K+].[#17-].[#53-]',
            'refined salt': '[Na+].[#17-].[#8-]-[Si](=[#8])-[#8-].[#8-]-[Si](=[#8])-[#8-].[Al+3].[#6](=[#8])(-[#8-])-[#8-].[Mg+2]',
            'sea salt': '[Na+].[K+].[Mg+2].[Ca+2].[#17-]',
            'smoked salt': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1.[Na+].[K+].[Mg+2].[Ca+2].[#17-]',
            'sodium iodide': '[Na+].[#53-]',
            'unrefined salt': '[Na+].[K+].[Mg+2].[Ca+2].[#17-].[Zn+2].[Cu+].[Fe+3].[#15]',
            'rock salt': '[Na+].[#17-]',
        }

        return smarts