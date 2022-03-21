#!/usr/bin/env python3
#
# GlobalChem - Bidendate Phosphine Ligands for Nickel
#
# ---------------------------------------------------

class NickelBidendatePhosphineLigands(object):

    def __init__(self):

        self.name = 'nickel_ligands'

    @staticmethod
    def get_smiles():

        smiles = {
            '1,1′-bis(diphenylphosphino)-methane': 'C1=CC=C(C=C1)P(CP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,1′-bis(dicyclohexylphosphino)-methane': 'C1CCC(CC1)P(CP(C2CCCCC2)C3CCCCC3)C4CCCCC4',
            '1,2-bis(diphenylphosphino)-ethane': 'C1=CC=C(C=C1)P(CCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,2-Bis(DI-tert-butylphosphino)ethane': 'CC(C)(C)P(CCP(C(C)(C)C)C(C)(C)C)C(C)(C)C',
            '1,2-bis(diisopropylphosphino)-ethane': 'CC(C)P(CCP(C(C)C)C(C)C)C(C)C',
            '1,2-bis(diethylphosphino)ethane': 'CCP(CC)CCP(CC)CC',
            '1,2-bis(dimethylphosphino)-ethane': 'CP(C)CCP(C)C',
            '1,2-bis(dicyclohexylphosphino)-ethane': 'C1CCC(CC1)P(CCP(C2CCCCC2)C3CCCCC3)C4CCCCC4',
            '1,3-Bis(diphenylphosphino)propane': 'C1=CC=C(C=C1)P(CCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,3-Bis(DI-I-propylphosphino)propane': 'CC(C)P(CCCP(C(C)C)C(C)C)C(C)C',
            '1,3-Bis(dimethylphosphino)propane': 'CP(C)CCCP(C)C',
            '1,3-bis(dicyclohexylphosphino)-propane': 'C1CCC(CC1)P(CCCP(C2CCCCC2)C3CCCCC3)C4CCCCC4',
            '1,4-bis(diphenylphosphino)-butane': 'C1=CC=C(C=C1)P(CCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,4-Bis(DI-I-propylphosphino)butane':'CC(C)P(CCCCP(C(C)C)C(C)C)C(C)C',
            '1,4-Bis(dicyclohexylphosphino)butane': 'C1CCC(CC1)P(CCCCP(C2CCCCC2)C3CCCCC3)C4CCCCC4',
            '1,4-bis(dicyclopentylphosphino)-butane': 'C1CC(CC1)P(CCCCP(C2CCCC2)C3CCCC3)C4CCCC4',
            '1,5-bis(diphenylphosphino)-pentane': 'C1=CC=C(C=C1)P(CCCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,6-Bis(diphenylphosphino)hexane': 'C1=CC=C(C=C1)P(CCCCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,8-bis(diphenylphosphino)-octane': 'C1=CC=C(C=C1)P(CCCCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,10-Bis(diphenylphosphino)decane': 'C1=CC=C(C=C1)P(CCCCCCCCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,1-Bis(diphenylphosphino)ethylene':'C=C(P(C1=CC=CC=C1)C2=CC=CC=C2)P(C3=CC=CC=C3)C4=CC=CC=C4',
            '1,2-bis(2,5-dimethylphospholan1-yl)ethane': 'CC1CCC(P1CCP2C(CCC2C)C)C.C1CC=CCCC=C1',
            '1,2-bis(2,5-diisopropylphospholan-1-yl)ethane': 'CC(C)C1CCC(P1CCP2C(CCC2C(C)C)C(C)C)C(C)C',
            '1,2-bis(2,5-diphenylphospholin1-yl)ethane': 'P(c1ccccc1)(c2ccccc2)CCP(c3ccccc3)c4ccccc4',
            'butane-2,3-diylbis(diphenylphosphine)': 'CC(C(C)P(C1=CC=CC=C1)C2=CC=CC=C2)P(C3=CC=CC=C3)C4=CC=CC=C4',
            '1,2-Bis(diphenylphosphino)propane':'CC(CP(C1=CC=CC=C1)C2=CC=CC=C2)P(C3=CC=CC=C3)C4=CC=CC=C4',
            '1,2-bis(phenyl(o-tolyl)-phosphino)ethane':'CC1=CC=CC=C1P(CCP(C2=CC=CC=C2)C3=CC=CC=C3C)C4=CC=CC=C4',
            '(Z)-1,2-Bis(diphenylphosphino)ethene':'C1=CC=C(C=C1)P(C=CP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4',
            '1,1′-ditert-butyl-2,2′-biphospholane': 'CC(C)OP(=O)(C(CC1=CC(=C(C(=C1)C(C)(C)C)O)C(C)(C)C)P(=O)(OC(C)C)OC(C)C)OC(C)C',
            '2,2′-ditert-butyl-2,2′,3,3′-tetrahydro-1H,1’H-1,1′-biisophosphindole':'CC(P1([CH]([CH]2P(CC3=C2C=CC=C3)(C(C)(C)C)=O)C4=C(C=CC=C4)C1)=O)(C)C',
            '2,3-bis( tert-butyl(methyl)-phosphino)quinoxaline': 'CC(C)(C)P(C)C1=NC2=CC=CC=C2N=C1P(C)C(C)(C)C',
            '3,4-bis(dicyclohexylphosphino)-thiophene': 'C1CCC(CC1)P(C2CCCCC2)C3=CSC=C3P(C4CCCCC4)C5CCCCC5',
            '(1R,4S,5R)-5,6-bis(diphenylphosphino)bicyclo[2.2.1]hept-2-ene': 'C1C2C=CC1C(C2P(C3=CC=CC=C3)C4=CC=CC=C4)P(C5=CC=CC=C5)C6=CC=CC=C6',
            '4,4′-di-tert-butyl-4,4′,5,5′-tetrahydro-3H,3′H-3,3′-bidinaphtho-[2,1-c:1′,2′-e]phosphepine': 'CC(C)(C)P1Cc2ccc3ccccc3c2-c4c(ccc5ccccc45)C1C6P(Cc7ccc8ccccc8c7-c9c6ccc%10ccccc9%10)C(C)(C)C',
            '1,2-Bis(diphenylphosphino)benzene':'C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3P(C4=CC=CC=C4)C5=CC=CC=C5',
            '1,2-bis(diethylphosphino)-benzene':'C1=CC=C(C=C1)P(CC)CCP(CC)CC',
            '1,2-bis(diisopropylphosphino)-benzene': 'C1=CC=C(C=C1)P(C(C)C)CCP(C(C)C)CC',
            '1,2-bis(dicyclohexylphosphino)-benzene': 'C1CCC(CC1)P(C2CCCCC2)C3=CC=CC=C3P(C4CCCCC4)C5CCCCC5',
            '1,2-Bis((R)-tert-butyl(methyl)phosphino)benzene': 'CC(C)(C)P(C)C1=CC=CC=C1P(C)C(C)(C)C',
            '1-(2-(2,5-dimethylphospholan-1-yl)phenyl)-2,3,5-trimethylphospholane':'CC1CCC(P1C2=CC=CC=C2P3C(CCC3C)C)C',
            '1-(2-(2,5-diethylphospholan-1-yl)phenyl)-2,3,5-trimethylphospholane':'CCC1CCC(P1C2=CC=CC=C2P3C(CCC3CC)CC)CC',
            '1-(2-(2,5-diisopropylphospholan1-yl)phenyl)-2,3,5-trimethylphospholane': 'CC(C)C1CCC(P1C2=CC=CC=C2P3C(CCC3C(C)C)C(C)C)C(C)C',
            'pentane-2,4-diylbis(diphenylphosphine)': 'CC(CC(C)P(C1=CC=CC=C1)C2=CC=CC=C2)P(C3=CC=CC=C3)C4=CC=CC=C4',
            '2,2′-bis(diphenylphosphino)-1,1′-binaphthalene': 'C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=C(C4=CC=CC=C4C=C3)C5=C(C=CC6=CC=CC=C65)P(C7=CC=CC=C7)C8=CC=CC=C8',
            '2,2′-bis(di(p-tol)phosphino)-1,1′-binaphthalene': 'Cc1ccc(cc1)P(c2ccc(C)cc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccc(C)cc7)c8ccc(C)cc8',
            '2,2′-bis(dicyclohexylphosphino)-1,1′-binaphthalene': 'C1CCC(CC1)P(C2CCCCC2)C3=C(C4=CC=CC=C4C=C3)C5=C(C=CC6=CC=CC=C65)P(C7CCCCC7)C8CCCCC8',
            '2,2′-bis(bis(3,5-dimethylphenyl)-phosphino)-1,1′-binaphthalene': 'Cc1cc(C)cc(c1)P(c2cc(C)cc(C)c2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7cc(C)cc(C)c7)c8cc(C)cc(C)c8',
            '5,5′-bis(diphenylphosphino)-,4′-bibenzo[d][1,3]dioxole':'C1OC2=C(O1)C(=C(C=C2)P(C3=CC=CC=C3)C4=CC=CC=C4)C5=C(C=CC6=C5OCO6)P(C7=CC=CC=C7)C8=CC=CC=C8',
            '6,6′-bis(diphenylphosphino)-2,2′,3,3′-tetrahydro-5,5′-bibenzo-[b][1,4]dioxine': 'C1COC2=C(O1)C=CC(=C2C3=C(C=CC4=C3OCCO4)P(C5=CC=CC=C5)C6=CC=CC=C6)P(C7=CC=CC=C7)C8=CC=CC=C8',
            '5,5′-bis(bis(3,5-dimethylphenyl)-phosphino)-4,4′-bibenzo[d]-[1,3]dioxole':'CC1=CC(=CC(=C1)P(C2=C(C3=C(C=C2)OCO3)C4=C(C=CC5=C4OCO5)P(C6=CC(=CC(=C6)C)C)C7=CC(=CC(=C7)C)C)C8=CC(=CC(=C8)C)C)C',
            '5,5′-bis(bis(3,5-ditert-butyl-4-methoxyphenyl)phosphino)-4,4′-bibenzo[d][1,3]dioxole': 'CC(C)(C)C1=CC(=CC(=C1OC)C(C)(C)C)P(C2=C(C3=C(C=C2)OCO3)C4=C(C=CC5=C4OCO5)P(C6=CC(=C(C(=C6)C(C)(C)C)OC)C(C)(C)C)C7=CC(=C(C(=C7)C(C)(C)C)OC)C(C)(C)C)C8=CC(=C(C(=C8)C(C)(C)C)OC)C(C)(C)C',
            '2,2′-bis(diphenylphosphino)-1,1′-biphenyl': 'C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3C4=CC=CC=C4P(C5=CC=CC=C5)C6=CC=CC=C6',
            '(6,6′-dimethyl-[1,1′-biphenyl]-2,2′-diyl)bis(diphenylphosphine)':'CC1=C(C(=CC=C1)P(C2=CC=CC=C2)C3=CC=CC=C3)C4=C(C=CC=C4P(C5=CC=CC=C5)C6=CC=CC=C6)C',
            '(6,6′-dimethoxy-[1,1′-biphenyl]-2,2′-diyl)bis(diphenylphosphane': 'COc1cccc(P(c2ccccc2)c3ccccc3)c1-c4c(OC)cccc4P(c5ccccc5)c6ccccc6',
            '(6,6′-dimethoxy-[1,1′-biphenyl]-2,2′-diyl)bis(bis(3,5-di-tert-butyl4-methoxyphenyl)phosphine)': 'COc1cccc(P(c2cc(c(OC)c(c2)C(C)(C)C)C(C)(C)C)c3cc(c(OC)c(c3)C(C)(C)C)C(C)(C)C)c1-c4c(OC)cccc4P(c5cc(c(OC)c(c5)C(C)(C)C)C(C)(C)C)c6cc(c(OC)c(c6)C(C)(C)C)C(C)(C)C',
            '(6,6′-dimethoxy-[1,1′-biphenyl]-2,2′-diyl)bis(bis(3,5-di-tertbutylphenyl)phosphine)': 'COc1cccc(P(c2cc(cc(c2)C(C)(C)C)C(C)(C)C)c3cc(cc(c3)C(C)(C)C)C(C)(C)C)c1-c4c(OC)cccc4P(c5cc(cc(c5)C(C)(C)C)C(C)(C)C)c6cc(cc(c6)C(C)(C)C)C(C)(C)C',

        }

        missing_entries = {
            '1,2-bis(bis(6,6,6,6,6-pentafluoro6l8-hexa-1,3,5-triyn-1-yl)-phosphino)ethane': '',
            '1,9-bis(diisopropylphosphino)-nonane': '',
            'benzo[b]thiophene-2,3-diylbis-(dicyclohexylphosphane)': '',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {

        }

        return smarts