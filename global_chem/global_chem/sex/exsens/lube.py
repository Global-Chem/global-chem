#!/usr/bin/env python3
#
# GlobalChem - Common Lube Ingredients
#
# ------------------------------------

class Lube(object):

    def __init__(self):

        self.name = 'lube'

    @staticmethod
    def get_smiles():

        smiles = {
            'water': 'O',
            'glycerin': 'C(C(CO)O)O',
            'dimethicone': 'C[Si](C)(C)O[Si](C)(C)C',
            'vinyl dimethicone': 'C[Si](C)(C)O[Si](C)(C)O[Si](C)(C=C)O[Si](C)(C)C',
            'cyclomethicone': 'C[Si]1(O[Si](O[Si](O[Si](O[Si](O1)(C)C)(C)C)(C)C)(C)C)C',
            'phenyl trimethicone': 'C[Si](C)(C)O[Si](C1=CC=CC=C1)(O[Si](C)(C)C)O[Si](C)(C)C',
            'cyclopentasiloxane': 'C[Si](O[Si]1(C)C)(O[Si](C)(O[Si](C)(C)O[Si](C)(C)O1)C)C',
            'polyethylene glycol': 'COCCCCO',
            'propylene glycol': 'CC(CO)O',
            'propanediol': 'CCC(O)O',
            'polyoxyethylene': 'CCCCCCCCC=CCCCCCCCC(=O)OCCOCC(C1C(CC(O1)OCCO)OCCO)OCCO',
            'optifio H370VF': '',
            'carboxymethylcellulose': 'CC(=O)[O-].C(C(C(C(C(C=O)O)O)O)O)O',
            'hydroxyethylcellulose': 'CC(COCC1C(C(C(C(O1)OC2C(OC(C(C2OCC(C)O)OCC(C)O)OCC(C)O)COCC(C)O)OCC(C)O)OCC(C)O)OCC(C)O)O',
            'xanthan gum': 'CC(=O)OCC1C(C(C(C(O1)OC2C(C(OC(C2OC3C(C(C(C(O3)CO)OP)O)O)CO)CP)O)OC4C(C(C(C(O4)C(=O)O)OC5C(C(C6C(O5)COC(O6)(C)C(=O)O)O)O)O)O)O)O',
            'agar': 'CC1C(C2C(C(O1)CO2)OC3C(C(C(C(O3)CO)O)OC)O)O',
            'pectin': 'C1(C(C(OC(C1O)O)C(=O)O)O)O',
            'maltodextrin': 'OCC(O)C(O)C(O)C(O)C=O',
            'sodium benzoate': 'C1=CC=C(C=C1)C(=O)[O-].[Na+]',
            'potassium sorbate': 'CC=CC=CC(=O)[O-].[K+]',
            'chlorhexidine': 'C1=CC(=CC=C1NC(=NC(=NCCCCCCN=C(N)N=C(N)NC2=CC=C(C=C2)Cl)N)N)Cl',
            'phenoxyethanol': 'C1=CC=C(C=C1)OCCO',
            'pethylhexylglycerin': 'OC(COCC(CCCC)CC)CO',
            'pentylene glycol': 'CCCC(CO)O',
            'methylparaben': 'COC(=O)C1=CC=C(C=C1)O',
            'propylparaben': 'CCCOC(=O)C1=CC=C(C=C1)O',
            'butylparaben': 'CCCCOC(=O)C1=CC=C(C=C1)O',
            'phthalates': 'COC(=O)C1=CC=CC=C1C(=O)OC',
            'chlorhexidine digluconate': 'C1=CC(=CC=C1NC(=NC(=NCCCCCCN=C(N)N=C(N)NC2=CC=C(C=C2)Cl)N)N)Cl.C(C(C(C(C(C(=O)O)O)O)O)O)O.C(C(C(C(C(C(=O)O)O)O)O)O)O',
            'diazolidinyl urea': 'C(NC(=O)N(CO)C1C(=O)N(C(=O)N1CO)CO)O',
            'sucralose': 'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CCl)O)O)CCl)O)O)Cl)O',
            'sodium saccharine': 'C1=CC=C2C(=C1)C(=O)NS2(=O)=O.[Na+',
            'rebaudioside a': 'CC12CCCC(C1CCC34C2CCC(C3)(C(=C)C4)OC5C(C(C(C(O5)CO)O)OC6C(C(C(C(O6)CO)O)O)O)OC7C(C(C(C(O7)CO)O)O)O)(C)C(=O)OC8C(C(C(C(O8)CO)O)O)O',
            'tocopheryl acetate': 'CC1=C(C(=C(C2=C1OC(CC2)(C)CCCC(C)CCCC(C)CCCC(C)C)C)OC(=O)C)C',
            'sodium hydroxide': 'NaO',
            'gluconolactone': 'C(C1C(C(C(C(=O)O1)O)O)O)O',
            'sodium hyaluronate': 'CC(=O)NC1CC(C(OC1OC2C(C(C(OC2C(=O)[O-])O)O)O)CO)O.[Na+]',
            'castor oil': 'CCCCCCC(CC=CCCCCCCCC(=O)OCC(COC(=O)CCCCCCCC=CCC(CCCCCC)O)OC(=O)CCCCCCCC=CCC(CCCCCC)O)O',
            'benzocaine': 'CCOC(=O)C1=CC=C(C=C1)N',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'water': '[#8]',
            'glycerin': '[#6](-[#6](-[#6]-[#8])-[#8])-[#8]',
            'dimethicone': '[#6]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#6]',
            'vinyl dimethicone': '[#6]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6])(-[#6]=[#6])-[#8]-[Si](-[#6])(-[#6])-[#6]',
            'cyclomethicone': '[#6]-[Si]1(-[#8]-[Si](-[#8]-[Si](-[#8]-[Si](-[#8]-[Si](-[#8]-1)(-[#6])-[#6])(-[#6])-[#6])(-[#6])-[#6])(-[#6])-[#6])-[#6]',
            'phenyl trimethicone': '[#6]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#8]-[Si](-[#6])(-[#6])-[#6])-[#8]-[Si](-[#6])(-[#6])-[#6]',
            'cyclopentasiloxane': '[#6]-[Si]1(-[#8]-[Si](-[#6])(-[#6])-[#8]-[Si](-[#8]-[Si](-[#8]-[Si](-[#8]-1)(-[#6])-[#6])(-[#6])-[#6])(-[#6])-[#6])-[#6]',
            'polyethylene glycol': '[#6]-[#8]-[#6]-[#6]-[#6]-[#6]-[#8]',
            'propylene glycol': '[#6]-[#6](-[#6]-[#8])-[#8]',
            'propanediol': '[#6]-[#6]-[#6](-[#8])-[#8]',
            'polyoxyethylene': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]=[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6]-[#8]-[#6]-[#6](-[#6]1-[#6](-[#6]-[#6](-[#8]-1)-[#8]-[#6]-[#6]-[#8])-[#8]-[#6]-[#6]-[#8])-[#8]-[#6]-[#6]-[#8]',
            'optifio H370VF': '',
            'carboxymethylcellulose': '[#6]-[#6](=[#8])-[#8-].[#6](-[#6](-[#6](-[#6](-[#6](-[#6]=[#8])-[#8])-[#8])-[#8])-[#8])-[#8]',
            'hydroxyethylcellulose': '[#6]-[#6](-[#6]-[#8]-[#6]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1-[#6](-[#8]-[#6](-[#6](-[#6]-1-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#6]-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]-[#6]-[#6](-[#6])-[#8])-[#8]',
            'xanthan gum': '[#6]-[#6](=[#8])-[#8]-[#6]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1-[#6](-[#6](-[#8]-[#6](-[#6]-1-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8]-[#15])-[#8])-[#8])-[#6]-[#8])-[#6]-[#15])-[#8])-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6](=[#8])-[#8])-[#8]-[#6]1-[#6](-[#6](-[#6]2-[#6](-[#8]-1)-[#6]-[#8]-[#6](-[#8]-2)(-[#6])-[#6](=[#8])-[#8])-[#8])-[#8])-[#8])-[#8])-[#8])-[#8]',
            'agar': '[#6]-[#6]1-[#6](-[#6]2-[#6](-[#6](-[#8]-1)-[#6]-[#8]-2)-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8]-[#6])-[#8])-[#8]',
            'pectin': '[#6]1(-[#6](-[#6](-[#8]-[#6](-[#6]-1-[#8])-[#8])-[#6](=[#8])-[#8])-[#8])-[#8]',
            'maltodextrin': '[#8]-[#6]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]=[#8]',
            'sodium benzoate': '[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](=[#8])-[#8-].[Na+]',
            'potassium sorbate': '[#6]-[#6]=[#6]-[#6]=[#6]-[#6](=[#8])-[#8-].[K+]',
            'chlorhexidine': '[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1-[#7]-[#6](=[#7]-[#6](=[#7]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#7]=[#6](-[#7])-[#7]=[#6](-[#7])-[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#17])-[#7])-[#7])-[#17]',
            'phenoxyethanol': '[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8]-[#6]-[#6]-[#8]',
            'pethylhexylglycerin': '[#8]-[#6](-[#6]-[#8]-[#6]-[#6](-[#6]-[#6]-[#6]-[#6])-[#6]-[#6])-[#6]-[#8]',
            'pentylene glycol': '[#6]-[#6]-[#6]-[#6](-[#6]-[#8])-[#8]',
            'methylparaben': '[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8]',
            'propylparaben': '[#6]-[#6]-[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8]',
            'butylparaben': '[#6]-[#6]-[#6]-[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8]',
            'phthalates': '[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](=[#8])-[#8]-[#6]',
            'chlorhexidine digluconate': '[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1-[#7]-[#6](=[#7]-[#6](=[#7]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#7]=[#6](-[#7])-[#7]=[#6](-[#7])-[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#17])-[#7])-[#7])-[#17].[#6](-[#6](-[#6](-[#6](-[#6](-[#6](=[#8])-[#8])-[#8])-[#8])-[#8])-[#8])-[#8].[#6](-[#6](-[#6](-[#6](-[#6](-[#6](=[#8])-[#8])-[#8])-[#8])-[#8])-[#8])-[#8]',
            'diazolidinyl urea': '[#6](-[#7]-[#6](=[#8])-[#7](-[#6]-[#8])-[#6]1-[#6](=[#8])-[#7](-[#6](=[#8])-[#7]-1-[#6]-[#8])-[#6]-[#8])-[#8]',
            'sucralose': '[#6](-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#8]-[#6]1(-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#17])-[#8])-[#8])-[#6]-[#17])-[#8])-[#8])-[#17])-[#8]',
            'sodium saccharine': '[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1)-[#6](=[#8])-[#7]-[#16]-2(=[#8])=[#8].[Na+]',
            'rebaudioside a': '[#6]-[#6]12-[#6]-[#6]-[#6]-[#6](-[#6]-1-[#6]-[#6]-[#6]13-[#6]-2-[#6]-[#6]-[#6](-[#6]-1)(-[#6](=[#6])-[#6]-3)-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8])-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8])(-[#6])-[#6](=[#8])-[#8]-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8]',
            'tocopheryl acetate': '[#6]-[#6]1:[#6](:[#6](:[#6](:[#6]2:[#6]:1-[#8]-[#6](-[#6]-[#6]-2)(-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6]-[#6](-[#6])-[#6])-[#6])-[#8]-[#6](=[#8])-[#6])-[#6]',
            'sodium hydroxide': '[Na+].[#8]',
            'gluconolactone': '[#6](-[#6]1-[#6](-[#6](-[#6](-[#6](=[#8])-[#8]-1)-[#8])-[#8])-[#8])-[#8]',
            'sodium hyaluronate': '[#6]-[#6](=[#8])-[#7]-[#6]1-[#6]-[#6](-[#6](-[#8]-[#6]-1-[#8]-[#6]1-[#6](-[#6](-[#6](-[#8]-[#6]-1-[#6](=[#8])-[#8-])-[#8])-[#8])-[#8])-[#6]-[#8])-[#8].[Na+]',
            'castor oil': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6]-[#6]=[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#8]-[#6]-[#6](-[#6]-[#8]-[#6](=[#8])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]=[#6]-[#6]-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])-[#8])-[#8]-[#6](=[#8])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]=[#6]-[#6]-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])-[#8])-[#8]',
            'benzocaine': '[#6]-[#6]-[#8]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#7]',
        }

        return smarts