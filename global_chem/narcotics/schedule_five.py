#!/usr/bin/env python3
#
# GlobalChem - Schedule Five
#
# -----------------------------------

class ScheduleFive(object):

    def __init__(self):

        pass

    def get_schedule_five_smiles(self):

        functional_groups = {
            'eluxadoline': 'CC1=CC(=CC(=C1CC(C(=O)N(CC2=CC(=C(C=C2)OC)C(=O)O)C(C)C3=NC=C(N3)C4=CC=CC=C4)N)C)C(=O)N',
            'pyrovalerone': 'CCCC(C(=O)C1=CC=C(C=C1)C)N2CCCC2',
            'brivaracetam': 'CCCC1CC(=O)N(C1)C(CC)C(=O)N',
            'cenobamate': 'C1=CC=C(C(=C1)C(CN2N=CN=N2)OC(=O)N)Cl',
            'ezogabine': 'CCOC(=O)NC1=C(C=C(C=C1)NCC2=CC=C(C=C2)F)N',
            'lacosamide': 'CC(=O)NC(COC)C(=O)NCC1=CC=CC=C1',
            'lasmiditan': 'CN1CCC(CC1)C(=O)C2=NC(=CC=C2)NC(=O)C3=C(C=C(C=C3F)F)F',
            'pregabalin': 'CC(C)CC(CC(=O)O)CN',
        }

    def get_schedule_five_smarts(self):

        functional_groups = {
            'eluxadoline': '[#6]-[#6]1:[#6]:[#6](:[#6]:[#6](:[#6]:1-[#6]-[#6](-[#6](=[#8])-[#7](-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6]:[#6]:1)-[#8]-[#6])-[#6](=[#8])-[#8])-[#6](-[#6])-[#6]1:[#7]:[#6]:[#6](:[#7H]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7])-[#6])-[#6](=[#8])-[#7]',
            'pyrovalerone': '[#6]-[#6]-[#6]-[#6](-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6])-[#7]1-[#6]-[#6]-[#6]-[#6]-1',
            'brivaracetam': '[#6]-[#6]-[#6]-[#6]1-[#6]-[#6](=[#8])-[#7](-[#6]-1)-[#6](-[#6]-[#6])-[#6](=[#8])-[#7]',
            'cenobamate': '[#6]1:[#6]:[#6]:[#6](:[#6](:[#6]:1)-[#6](-[#6]-[#7]1:[#7]:[#6]:[#7]:[#7]:1)-[#8]-[#6](=[#8])-[#7])-[#17]',
            'ezogabine': '[#6]-[#6]-[#8]-[#6](=[#8])-[#7]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1)-[#7]-[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#9])-[#7]',
            'lacosamide': '[#6]-[#6](=[#8])-[#7]-[#6](-[#6]-[#8]-[#6])-[#6](=[#8])-[#7]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'lasmiditan': '[#6]-[#7]1-[#6]-[#6]-[#6](-[#6]-[#6]-1)-[#6](=[#8])-[#6]1:[#7]:[#6](:[#6]:[#6]:[#6]:1)-[#7]-[#6](=[#8])-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#9])-[#9])-[#9]',
            'pregabalin': '[#6]-[#6](-[#6])-[#6]-[#6](-[#6]-[#6](=[#8])-[#8])-[#6]-[#7]',
        }