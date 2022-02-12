#!/usr/bin/env python3
#
# GlobalChem - Warheads for Covalent Inhibitors
#
# ---------------------------------------------

class CommonWarheadsCovalentInhibitors(object):

    def __init__(self):

        pass

    @staticmethod
    def get_smiles():

        smiles = {
            'propiolamide': 'C#CC(N)=O',
            'fumarate ester': 'NC(/C=C/CC(OC)=O)=O',
            'allenamide': 'NC(C=C=C)=O',
            'propiolonitrile': 'C#CC#N',
            'propargylamide': 'C#CCC(N)=O',
            'arylsulfonyl bicyclobutane': 'O=S(C12CC1C2)(C3=CC=CC=C3)=O',
            'haloalkane': 'CBr',
            'alpha-halomethyl': 'CC(CCl)=O',
            'alpha-haloamide': 'NC(CCl)=O',
            'alpha-haloester': 'O=C(CCl)OC',
            'epoxide': 'C1CO1',
            'aziridine': 'N1CC1',
            'nitroalkane': 'CC[N+]([O-])=O',
            'acrylamide': 'C=CC(N)=O',
            'cyanoenone': 'O=C(C)C(C#N)=C',
            'aldehyde': 'O=C(C)[H]',
            'ketone': 'O=C(C)C',
            'nitrile': 'N#CC',
            'cyanamide': 'NC#N',
            'isothicyanate': '[N-]=C=S',
            'sulfone': 'CS=O',
            'sulfonyl fluoride': 'O=S(F)=O',
            'sulfonimidoyl fluoride': 'N=S(F)(F)=O',
            'aryl fluorosulfate': 'O=S(OCCCCC)(F)=O',
            'ester': 'CC(OC)=O',
            'sulfonamide': 'O=S(N)=O',
            '2-carbonyl arylboronic acid': 'O=C(C1=CC=CC=C1B(O)O)C',
            'n-methyl isoxazolium': 'C[N+]1=CC=CO1',
            'oxaziridine': 'O1NC1',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'propiolamide': '[#6]#[#6]-[#6](-[#7])=[#8]',
            'fumarate ester': '[#7]-[#6](/[#6]=[#6]/[#6]-[#6](-[#8]-[#6])=[#8])=[#8]',
            'allenamide': '[#7]-[#6](-[#6]=[#6]=[#6])=[#8]',
            'propiolonitrile': '[#6]#[#6]-[#6]#[#7]',
            'propargylamide': '[#6]#[#6]-[#6]-[#6](-[#7])=[#8]',
            'arylsulfonyl bicyclobutane': '[#8]=[#16](-[#6]12-[#6]-[#6]-1-[#6]-2)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8]',
            'haloalkane': '[#6]-[#35]',
            'alpha-halomethyl': '[#6]-[#6](-[#6]-[#17])=[#8]',
            'alpha-haloamide': '[#7]-[#6](-[#6]-[#17])=[#8]',
            'alpha-haloester': '[#8]=[#6](-[#6]-[#17])-[#8]-[#6]',
            'epoxide': '[#6]1-[#6]-[#8]-1',
            'aziridine': '[#7]1-[#6]-[#6]-1',
            'nitroalkane': '[#6]-[#6]-[#7+](-[#8-])=[#8]',
            'acrylamide': '[#6]=[#6]-[#6](-[#7])=[#8]',
            'cyanoenone': '[#8]=[#6](-[#6])-[#6](-[#6]#[#7])=[#6]',
            'aldehyde': '[#8]=[#6H]-[#6]',
            'ketone': '[#8]=[#6](-[#6])-[#6]',
            'nitrile': '[#7]#[#6]-[#6]',
            'cyanamide': '[#7]-[#6]#[#7]',
            'isothicyanate': '[#7-]=[#6]=[#16]',
            'sulfone': '[#6]-[#16]=[#8]',
            'sulfonyl fluoride': '[#8]=[#16](-[#9])=[#8]',
            'sulfonimidoyl fluoride': '[#7]=[#16](-[#9])(-[#9])=[#8]',
            'aryl fluorosulfate': '[#8]=[#16](-[#8]-[#6]-[#6]-[#6]-[#6]-[#6])(-[#9])=[#8]',
            'ester': '[#6]-[#6](-[#8]-[#6])=[#8]',
            'sulfonamide': '[#8]=[#16](-[#7])=[#8]',
            '2-carbonyl arylboronic acid': '[#8]=[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#5](-[#8])-[#8])-[#6]',
            'n-methyl isoxazolium': '[#6]-[#7+]1:[#6]:[#6]:[#6]:[#8]:1',
            'oxaziridine': '[#8]1-[#7]-[#6]-1',
        }

        return smarts
