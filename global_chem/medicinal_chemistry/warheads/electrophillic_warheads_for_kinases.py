#!/usr/bin/env python3
#
# GlobalChem - Electrophillic Warheads For Kinases
#
# ------------------------------------------------

class ElectrophilicWarheadsForKinases(object):

    def __init__(self):

        pass

    @staticmethod
    def get_smiles():

        smiles = {
            'methylacrylamide': 'CNC(C=C)=O',
            'methyl acrylate': 'COC(C=C)=O',
            'methyl propiolate' :'COC(C#C)=O',
            '2-cyanoacrylamide': 'N#CC(C(N)=O)=C',
            'n-methylmaleimide' :'CN1C(C=CC1=O)=O',
            'n-ethylmaleimide' :'O=C(C=CC1=O)N1CC',
            'crotonamide': 'C/C=C/C(N)=O',
            'ethyl crotonate': 'C/C=C/C(OCC)=O',
            'crotononitrile' :'C/C=C/C#N',
            'methyl methylpropiolate': 'CC#CC(OC)=O',
            'isothiocyanatomethane': 'CN=C=S',
            'isothiocyanatoethane': 'CCN=C=S',
            'prop-1-ene': 'CC=C',
            'prop-1-yne': 'CC#C',
            'acetonitrile': 'CC#N',
            'tert-butyl (Z)-2-ethylidenehydrazine-1-carboxylate' : 'C/C=N/\/NC(OC(C)(C)C)=O',
            'n-methylchloroacetamide': 'CNC(CCl)=O',
            'n-methyl-2-chloropropanamide': 'CNC(C(C)Cl)=O',
            'n-methyl-2-bromopropanamide': 'CNC(C(C)Br)=O',
            'bromoacetone': 'CC(CBr)=O',
            '2-methyloxirane': 'CC1OC1',
            'fluoromethane': 'CF',
            'methylsulfane': 'CS',
            'aldehyde': 'CC=O'
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'methylacrylamide': '[#6]-[#7]-[#6](-[#6]=[#6])=[#8]',
            'methyl acrylate': '[#6]-[#8]-[#6](-[#6]=[#6])=[#8]',
            'methyl propiolate': '[#6]-[#8]-[#6](-[#6]#[#6])=[#8]',
            '2-cyanoacrylamide': '[#7]#[#6]-[#6](-[#6](-[#7])=[#8])=[#6]',
            'n-methylmaleimide': '[#6]-[#7]1-[#6](-[#6]=[#6]-[#6]-1=[#8])=[#8]',
            'n-ethylmaleimide': '[#8]=[#6]1-[#6]=[#6]-[#6](=[#8])-[#7]-1-[#6]-[#6]',
            'crotonamide': '[#6]/[#6]=[#6]/[#6](-[#7])=[#8]',
            'ethyl crotonate': '[#6]/[#6]=[#6]/[#6](-[#8]-[#6]-[#6])=[#8]',
            'crotononitrile': '[#6]/[#6]=[#6]/[#6]#[#7]',
            'methyl methylpropiolate': '[#6]-[#6]#[#6]-[#6](-[#8]-[#6])=[#8]',
            'isothiocyanatomethane': '[#6]-[#7]=[#6]=[#16]',
            'isothiocyanatoethane': '[#6]-[#6]-[#7]=[#6]=[#16]',
            'prop-1-ene': '[#6]-[#6]=[#6]',
            'prop-1-yne': '[#6]-[#6]#[#6]',
            'acetonitrile': '[#6]-[#6]#[#7]',
            'tert-butyl (Z)-2-ethylidenehydrazine-1-carboxylate': '[#6]/[#6]=[#7]\[#7]-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])=[#8]',
            'n-methylchloroacetamide': '[#6]-[#7]-[#6](-[#6]-[#17])=[#8]',
            'n-methyl-2-chloropropanamide': '[#6]-[#7]-[#6](-[#6](-[#6])-[#17])=[#8]',
            'n-methyl-2-bromopropanamide': '[#6]-[#7]-[#6](-[#6](-[#6])-[#35])=[#8]',
            'bromoacetone': '[#6]-[#6](-[#6]-[#35])=[#8]',
            '2-methyloxirane': '[#6]-[#6]1-[#8]-[#6]-1',
            'fluoromethane': '[#6]-[#9]',
            'methylsulfane': '[#6]-[#16]',
            'aldehyde': '[#6]-[#6]=[#8]',
        }

        return smarts
