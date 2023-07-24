#!/usr/bin/env python3
#
# GlobalChem - Natural Fibers
#
# ---------------------------------------

class Natural_Fibers(object):

    def __init__(self):

        self.name = 'natural_fibers'

    @staticmethod
    def get_smiles():
       
        smiles = {
        'p-coumaric acid': 'C1=CC(=CC=C1C=CC(=O)O)O',
        'ferulic acid': 'COC1=C(C=CC(=C1)C=CC(=O)O)O',
        'syringyl alcohol': 'CC(=O)OCC1=CC(=C(C(=C1)OC)OC(=O)C)OC',
        'gallic acid': 'C1=C(C=C(C(=C1O)O)O)C(=O)O',
        'cellulose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
        'beta glucose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'starch': '', # What type of starch
        'alpha glucose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'hemicellulose': 'CC(C(C1=CC=CC=C1)O)NC.CC(C(C1=CC=CC=C1)O)NC',
        'arabinose': 'C1C(C(C(C(O1)O)O)O)O',
        'mannose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'glucose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'xylose': 'C1C(C(C(C(O1)O)O)O)O',        
        }

        return smiles