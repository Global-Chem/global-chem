#!/usr/bin/env python3
#
# GlobalChem - Chemicals From Biomass
#
# -----------------------------------

class ChemicalsFromBioMass(object):

    def __init__(self):

        self.name = 'chemicals_from_biomass'

    @staticmethod
    def get_smiles():

        '''

        Missing Entries:
             '8:2 chlorinated polyfluorinated ether sulfonic acid': '',

        '''

        smiles = {
            'synthesis gas': '[C-]#[O+].[HH]',
            'ethane': 'CC',
            'ethylene': 'C=C',
            'propylene': 'CC=C',
            '1,3-butadiene': 'C=CC=C',
            'butene': 'CCC=C',
            'butane': 'CCCC',
            'dicyclopentadiene': 'C1C=CC2C1C3CC2C=C3',
            'isoprene': 'CC(=C)C=C',
            'pentene': 'CCCC=C',
            'piperylene': 'CC=CC=C',
            'benzene': 'C1=CC=CC=C1',
            'toluene': 'CC1=CC=CC=C1',
            'p-xylene': 'CC1=CC=C(C=C1)C',
            'o-xylene': 'CC1=CC=CC=C1C',
            'm-xylene': 'CC1=CC(=CC=C1)C',
            'glycerol': 'C(C(CO)O)O',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'synthesis gas':'[#6-]#[#8+].[HH]',
            'ethane':'[#6]-[#6]',
            'ethylene':'[#6]=[#6]',
            'propylene':'[#6]-[#6]=[#6]',
            '1,3-butadiene':'[#6]=[#6]-[#6]=[#6]',
            'butene':'[#6]-[#6]-[#6]=[#6]',
            'butane':'[#6]-[#6]-[#6]-[#6]',
            'dicyclopentadiene':'[#6]1-[#6]=[#6]-[#6]2-[#6]-1-[#6]1-[#6]-[#6]-2-[#6]=[#6]-1',
            'isoprene':'[#6]-[#6](=[#6])-[#6]=[#6]',
            'pentene':'[#6]-[#6]-[#6]-[#6]=[#6]',
            'piperylene':'[#6]-[#6]=[#6]-[#6]=[#6]',
            'benzene':'[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'toluene':'[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'p-xylene':'[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]',
            'o-xylene':'[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]',
            'm-xylene':'[#6]-[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1)-[#6]',
            'glycerol':'[#6](-[#6](-[#6]-[#8])-[#8])-[#8]',
        }

        return smarts