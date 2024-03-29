#!/usr/bin/env python3
#
# GlobalChem - Schedule Two
#
# -----------------------------------

class ScheduleThree(object):

    def __init__(self):

        self.name = 'schedule_three'

    @staticmethod
    def get_smiles():

        smiles = {
            'norfentanyl': 'CCC(=O)N(C1CCNCC1)C2=CC=CC=C2',
            'benzphetamine': 'CC(CC1=CC=CC=C1)N(C)CC2=CC=CC=C2',
            'chlorphentermine': 'CN(C)CCC(C1=CC=C(C=C1)Cl)C2=CC=CC=N2',
            'clortermine': 'CC(C)(CC1=CC=CC=C1Cl)N',
            'phendimetrazine': 'CC1C(OCCN1C)C2=CC=CC=C2',
            'amobarbital': 'CCC1(C(=O)NC(=O)NC1=O)CCC(C)C',
            'secobarbital': 'CCCC(C)C1(C(=O)NC(=O)NC1=O)CC=C',
            'pentobarbital': 'CCCC(C)C1(C(=O)NC(=O)NC1=O)CC',
            'chlorhexadol': 'CC(CC(C)(C)O)OC(C(Cl)(Cl)Cl)O',
            'embutramide': 'CCC(CC)(CNC(=O)CCCO)C1=CC(=CC=C1)OC',
            'ketamine': 'CNC1(CCCCC1=O)C2=CC=CC=C2Cl',
            'lysergicacid': 'CN1CC(C=C2C1CC3=CNC4=CC=CC2=C34)C(=O)O',
            'lysergicacidamide': 'CN1CC(C=C2C1CC3=CNC4=CC=CC2=C34)C(=O)N',
            'methyprylon': 'CCC1(C(=O)C(CNC1=O)C)CC',
            'perampanel': 'C1=CC=C(C=C1)N2C=C(C=C(C2=O)C3=CC=CC=C3C#N)C4=CC=CC=N4',
            'sulfondiethylmethane': 'CCC(CC)(S(=O)(=O)CC)S(=O)(=O)CC',
            'sulfonethylmethane': 'CCC(C)(S(=O)(=O)CC)S(=O)(=O)CC',
            'sulfonmethane': 'O=S(=O)(C(C)(C)S(=O)(=O)CC)CC',
            'tiletamine': 'CCNC1(CCCCC1=O)C2=CC=CS2',
            'nalorphine_9400': 'C=CCN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',
            'flupyrazapon': 'CC1=NN(C2=C1C(=NCC(=O)N2C)C3=CC=CC=C3F)C',
            'buprenorphine': 'CC(C)(C)C(C)(C1CC23CCC1(C4C25CCN(C3CC6=C5C(=C(C=C6)O)O4)CC7CC7)OC)O',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'norfentanyl': '[#6]-[#6]-[#6](=[#8])-[#7](-[#6]1-[#6]-[#6]-[#7]-[#6]-[#6]-1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'benzphetamine': '[#6]-[#6](-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7](-[#6])-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'chlorphentermine': '[#6]-[#7](-[#6])-[#6]-[#6]-[#6](-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#17])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1',
            'clortermine': '[#6]-[#6](-[#6])(-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#17])-[#7]',
            'phendimetrazine': '[#6]-[#6]1-[#6](-[#8]-[#6]-[#6]-[#7]-1-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'amobarbital': '[#6]-[#6]-[#6]1(-[#6](=[#8])-[#7]-[#6](=[#8])-[#7]-[#6]-1=[#8])-[#6]-[#6]-[#6](-[#6])-[#6]',
            'secobarbital': '[#6]-[#6]-[#6]-[#6](-[#6])-[#6]1(-[#6](=[#8])-[#7]-[#6](=[#8])-[#7]-[#6]-1=[#8])-[#6]-[#6]=[#6]',
            'pentobarbital': '[#6]-[#6]-[#6]-[#6](-[#6])-[#6]1(-[#6](=[#8])-[#7]-[#6](=[#8])-[#7]-[#6]-1=[#8])-[#6]-[#6]',
            'chlorhexadol': '[#6]-[#6](-[#6]-[#6](-[#6])(-[#6])-[#8])-[#8]-[#6](-[#6](-[#17])(-[#17])-[#17])-[#8]',
            'embutramide': '[#6]-[#6]-[#6](-[#6]-[#6])(-[#6]-[#7]-[#6](=[#8])-[#6]-[#6]-[#6]-[#8])-[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1)-[#8]-[#6]',
            'ketamine': '[#6]-[#7]-[#6]1(-[#6]-[#6]-[#6]-[#6]-[#6]-1=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#17]',
            'lysergicacid': '[#6]-[#7]1-[#6]-[#6](-[#6]=[#6]2-[#6]-1-[#6]-[#6]1:[#6]:[#7H]:[#6]3:[#6]:[#6]:[#6]:[#6]-2:[#6]:1:3)-[#6](=[#8])-[#8]',
            'lysergicacidamide': '[#6]-[#7]1-[#6]-[#6](-[#6]=[#6]2-[#6]-1-[#6]-[#6]1:[#6]:[#7H]:[#6]3:[#6]:[#6]:[#6]:[#6]-2:[#6]:1:3)-[#6](=[#8])-[#7]',
            'methyprylon': '[#6]-[#6]-[#6]1(-[#6](=[#8])-[#6](-[#6]-[#7]-[#6]-1=[#8])-[#6])-[#6]-[#6]',
            'perampanel': '[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#7]1:[#6]:[#6](:[#6]:[#6](:[#6]:1=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]#[#7])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1',
            'sulfondiethylmethane': '[#6]-[#6]-[#6](-[#6]-[#6])(-[#16](=[#8])(=[#8])-[#6]-[#6])-[#16](=[#8])(=[#8])-[#6]-[#6]',
            'sulfonethylmethane': '[#6]-[#6]-[#6](-[#6])(-[#16](=[#8])(=[#8])-[#6]-[#6])-[#16](=[#8])(=[#8])-[#6]-[#6]',
            'sulfonmethane': '[#8]=[#16](=[#8])(-[#6](-[#6])(-[#6])-[#16](=[#8])(=[#8])-[#6]-[#6])-[#6]-[#6]',
            'tiletamine': '[#6]-[#6]-[#7]-[#6]1(-[#6]-[#6]-[#6]-[#6]-[#6]-1=[#8])-[#6]1:[#6]:[#6]:[#6]:[#16]:1',
            'nalorphine_9400': '[#6]=[#6]-[#6]-[#7]1-[#6]-[#6]-[#6]23-[#6]4-[#6]-1-[#6]-[#6]1:[#6]-2:[#6](:[#6](:[#6]:[#6]:1)-[#8])-[#8]-[#6]-3-[#6](-[#6]=[#6]-4)-[#8]',
            'flupyrazapon': '[#6]-[#6]1:[#7]:[#7](:[#6]2:[#6]:1-[#6](=[#7]-[#6]-[#6](=[#8])-[#7]-2-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#9])-[#6]',
            'buprenorphine': '[#6]-[#6](-[#6])(-[#6])-[#6](-[#6])(-[#6]1-[#6]-[#6]23-[#6]-[#6]-[#6]-1(-[#6]1-[#6]-24-[#6]-[#6]-[#7](-[#6]-3-[#6]-[#6]2:[#6]-4:[#6](:[#6](:[#6]:[#6]:2)-[#8])-[#8]-1)-[#6]-[#6]1-[#6]-[#6]-1)-[#8]-[#6])-[#8]',
        }

        return smarts