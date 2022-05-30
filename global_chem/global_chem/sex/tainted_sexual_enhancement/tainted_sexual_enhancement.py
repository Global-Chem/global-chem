#!/usr/bin/env python3
#
# GlobalChem - Common Ingredients for Tainted Sexual Enhancement
#
# --------------------------------------------------------------

class TaintedSexualEnhancements(object):

    def __init__(self):

        self.name = 'tainted_sexual_enhancements'

    @staticmethod
    def get_smiles():

        smiles = {
            'sildenafil': 'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
            'tadalafil,': 'CN1CC(=O)N2C(C1=O)CC3=C(C2C4=CC5=C(C=C4)OCO5)NC6=CC=CC=C36',
            'nitroglycerin': 'C(C(CO[N+](=O)[O-])O[N+](=O)[O-])O[N+](=O)[O-]',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'sildenafil': '[#6]-[#6]-[#6]-[#6]1:[#7]:[#7](:[#6]2:[#6]:1:[#7]:[#6](:[#7H]:[#6]:2=[#8])-[#6]1:[#6](:[#6]:[#6]:[#6](:[#6]:1)-[#16](=[#8])(=[#8])-[#7]1-[#6]-[#6]-[#7](-[#6]-[#6]-1)-[#6])-[#8]-[#6]-[#6])-[#6]',
            'tadalafil,': '[#6]-[#7]1-[#6]-[#6](=[#8])-[#7]2-[#6](-[#6]-1=[#8])-[#6]-[#6]1:[#6](-[#6]-2-[#6]2:[#6]:[#6]3:[#6](:[#6]:[#6]:2)-[#8]-[#6]-[#8]-3):[#7H]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:1:2',
            'nitroglycerin': '[#6](-[#6](-[#6]-[#8]-[#7+](=[#8])-[#8-])-[#8]-[#7+](=[#8])-[#8-])-[#8]-[#7+](=[#8])-[#8-]',
        }

        return smarts
