#!/usr/bin/env python3
#
# GlobalChem - Cannabinoids
#
# -----------------------------------

class Cannabinoids(object):

    def __init__(self):

        self.name = 'cannabinoids'

    @staticmethod
    def get_smiles():

        smiles = {
            'cannabigerolid acid': r'C/C(C)=C\CC/C(C)=C/CC1=C(O)C=C(C2=CC=CC=C2)C(C(O)=O)=C1O',
            'cannabigerolid acid monomethylether': r'C/C(C)=C\CC/C(C)=C/CC1=C(OC)C=C(C2=CC=CC=C2)C(C(O)=O)=C1O',
            'cannabigerol': r'C/C(C)=C\CC/C(C)=C/CC1=C(O)C=C(C2=CC=CC=C2)C=C1O',
            'cannabigerol monomethylether': r'C/C(C)=C\CC/C(C)=C/CC1=C(OC)C=C(C2=CC=CC=C2)C=C1O',
            'cannabigerovarinic acid': 'CCCC1=CC(=C(C(=C1C(=O)O)O)CC=C(C)CCC=C(C)C)O',
            'cannabigerovarin': 'CCCC1=CC(=C(C(=C1)O)CC=C(C)CCC=C(C)C)O',
            'cannabichromenic acid': 'CCCCCC1=CC2=C(C=CC(O2)(C)CCC=C(C)C)C(=C1C(=O)O)O',
            'cannabichromene': 'CCCCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
            'cannabichromevarinic acid': 'CCCC1=CC2=C(C=CC(O2)(C)CCC=C(C)C)C(=C1C(=O)O)O',
            'cannabichromevarin': 'CCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
            'cannabidiolic acid': 'CCCCCC1=CC(=C(C(=C1C(=O)O)O)C2C=C(CCC2C(=C)C)C)O',
            'cannabidiol': 'CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O',
            'cannabidiol monomethylether': 'CCCCCC1=CC(=C(C(=C1)OC)C2C=C(CCC2C(=C)C)C)O',
            'cannabidiol-C4': 'CC(CC[C@H]1C(C)=C)=CC1C2=CC=C(CCCC)C=C2O',
            'cannabidivarinic acid': 'CCCC1=CC(=C(C(=C1C(=O)O)O)C2C=C(CCC2C(=C)C)C)O',
            'cannabidivarin': 'CCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O',
            'cannabidiorcol': 'CC(CC[C@H]1C(C)=C)=CC1C2=CC=C(C)C=C2O',
            'delta-9-(trans)-tetrahydrocannabinolic acid A': 'CC(CC1)=CC(C1C(C)(C)O2)C3=C2C=C(C4=CC=CC=C4)C(C(O)=O)=C3O',
            'delta-9-(trans)-tetrahydrocannabinolic acid B': 'CC(CC1)=CC(C1C(C)(C)O2)C3=C2C(C(O)=O)=C(C4=CC=CC=C4)C=C3O',
            'delta-9-(trans)-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
            'delta-9-(trans)-tetrahydrocannabinolic acid C4': '',
            'delta-9-(trans)-tetrahydrocannabinol-C4': '',
            'delta-9-(trans)-tetrahydrocannabivarinic acid': '',
            'delta-9-(trans)-tetrahydrocannabiorcolic acid': '',
            'delta-9-(trans)-tetrahydrocannabiorcol': '',
            'delta-8-(trans)-tetrahydrocannabinolic acid': '',
            'delta-8-(trans)-tetrahydrocannabinol': '',
            'cannabicyclolic acid': '',
            'cannabicyclol': '',
            'cannabicyclovarin': '',
            'cannabielsoic acid A': '',
            'cannabielsoic acid B': '',
            'cannabielsoin': '',
            'cannabinolic acid': '',
            'cannabinol': '',
            'cannabinol methylether': '',
            'cannabinol C4': '',
            'cannabivarin': '',
            'cannabiorcol': '',
            '(-)-Cannabitriol': '',
            '(+)-Cannabitriol': '',
            '(=)-9,10-dihvdroxy-delta-6a-(10a)-tetrahydrocannabinol': '',
            '(-)-10-ethoxv-9-hydroxv-6a-(10a)-tetrahydrocannabinol': '',
            '(=)-8,9-dihvdroxy-6a-(10a)-tetrahvdrocannabinol': '',
            'cannabidiolic acid tetrahydrocannabitriol ester': '',
            'dehydrocannabifuran': '',
            'cannabifuran': '',
            'cannabichromanon': '',
            'cannabicitran': '',
            '10-oxo-delta-6a-(10a)-tetrahydrocannabinol': '',
            'delta-9-(6a,10a,cis)-tetrahydrocannabinol': '',
            '3,4,5,6-Tetrahydro-7-hydroxy-a,a-2-trimethyl9-n-propyl-2,6-methano-2H-l-benzoxocin-5-methanol': '',
            'cannabiripsol': '',
            '6a,7,1Oa-trihydroxy-delta-9-tetrahydrocannabinol': '',
            'choline': '',
            'trigonelline': '',
            'muscarine': '',
            'L-(+)-isoleucine betaine': '',
            'neurine': '',

        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {

        }

        return smarts