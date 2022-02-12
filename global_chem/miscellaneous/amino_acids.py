#!/usr/bin/env python3
#
# GlobalChem - Amino Acids
#
# -----------------------------------

class AminoAcids(object):

    def __init__(self):

        pass

    @staticmethod
    def get_smiles():

        smiles = {
            "alanine": "C",
            "arginine": "CCCCNC(N)=N",
            "asparagine": "CCC(N)=O",
            "aspartic acid": "CC(O)=O",
            "cysteine": "CS",
            "glutamic acid": "CCC(O)=O",
            "glutamine": "CCC(N)=O",
            "glycine": "[H]",
            "histidine": "CC1=CNC=N1",
            "isoleucine": "C(CC)([H])C",
            "leucine": "CC(C)C",
            "lysine": "CCCCN",
            "methionine": "CCSC",
            "phenylalanine": "CC1=CC=CC=C1",
            "proline": "C2CCCN2",
            "serine": "CO",
            "threonine": "C(C)([H])O",
            "tryptophan": "CCC1=CNC2=C1C=CC=C2",
            "tyrosine": "CC1=CC=C(O)C=C1",
            "valine": "C(C)C"
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'alanine':'[#6]',
            'arginine':'[#6]-[#6]-[#6]-[#6]-[#7]-[#6](-[#7])=[#7]',
            'asparagine':'[#6]-[#6]-[#6](-[#7])=[#8]',
            'aspartic acid':'[#6]-[#6](-[#8])=[#8]',
            'cysteine':'[#6]-[#16]',
            'glutamic acid':'[#6]-[#6]-[#6](-[#8])=[#8]',
            'glutamine':'[#6]-[#6]-[#6](-[#7])=[#8]',
            'glycine':'[H]',
            'histidine':'[#6]-[#6]1:[#6]:[#7H]:[#6]:[#7]:1',
            'isoleucine':'[#6H](-[#6]-[#6])-[#6]',
            'leucine':'[#6]-[#6](-[#6])-[#6]',
            'lysine':'[#6]-[#6]-[#6]-[#6]-[#7]',
            'methionine':'[#6]-[#6]-[#16]-[#6]',
            'phenylalanine':'[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'proline':'[#6]1-[#6]-[#6]-[#6]-[#7]-1',
            'serine':'[#6]-[#8]',
            'threonine':'[#6H](-[#6])-[#8]',
            'tryptophan':'[#6]-[#6]-[#6]1:[#6]:[#7H]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
            'tyrosine':'[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]):[#6]:[#6]:1',
            'valine':'[#6](-[#6])-[#6]',
        }

        return smarts
