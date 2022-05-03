#!/usr/bin/env python3
#
# GlobalChem - I Want to Live Longer
#
# -----------------------------------

class HowToLiveLonger(object):

    def __init__(self):

        self.name = 'i_want_to_live_longer'


    @staticmethod
    def get_smiles():

        smiles = {
            'metformin': 'CN(C)C(=N)N=C(N)N',
            'glucosamine': 'NC1C(OC(C(O)C1O)CO)O',
            'chondroitin_sulfate': 'CC(=O)NC1C(C(C(OC1O)OS(=O)(=O)O)O)OC2C(C(C(C(O2)C(=O)O)O)O)O',
            'spermidine': 'C(CCNCCCN)CN',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'metformin': '[#6]-[#7](-[#6])-[#6](=[#7])-[#7]=[#6](-[#7])-[#7]',
            'glucosamine': '[#7]-[#6]1-[#6](-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8])-[#6]-[#8])-[#8]',
            'chondroitin_sulfate': '[#6]-[#6](=[#8])-[#7]-[#6]1-[#6](-[#6](-[#6](-[#8]-[#6]-1-[#8])-[#8]-[#16](=[#8])(=[#8])-[#8])-[#8])-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6](=[#8])-[#8])-[#8])-[#8])-[#8]',
            'spermidine': '[#6](-[#6]-[#6]-[#7]-[#6]-[#6]-[#6]-[#7])-[#6]-[#7]',
        }

        return smarts