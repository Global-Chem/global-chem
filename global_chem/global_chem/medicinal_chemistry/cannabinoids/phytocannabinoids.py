#!/usr/bin/env python3
#
# GlobalChem - Cannabinoids
#
# -----------------------------------

class PhytoCannabinoids(object):

    def __init__(self):

        self.name = 'phytocannabinoids'

    @staticmethod
    def get_smiles():

        smiles = {
            'cannabinerol': r'CCCCCC1=CC(=C(C(=C1)O)C/C=C(/C)\CCC=C(C)C)O',
            'cannabinerolic acid': r'CCCCCC1=CC(=C(C(=C1C(=O)O)O)C/C=C(/C)\CCC=C(C)C)O',
            'carmagerol': r'CCCCCC1=CC(=C(C(=C1)O)C/C=C(\C)/CCC(C(C)(C)O)O)O',
            'rac-6-epoxycannabigerol': r'CCCCCC1=CC(=C(C(=C1)O)CC(O2)C2(CC\C=C(C)/C)C)O',
            'rac-6-epoxycannabigerolic acid': r'CCCCCC1=CC(=C(C(=C1C(=O)O)O)CC(O2)C2(CC\C=C(C)/C)C)O',
            'rac-6-epoxycannabinerol': r'CCCCCC1=CC(=C(C(=C1)O)CC(O2)C2(CC/C=C(C)\C)C)O',
            'rac-6-epoxycannabinerolic acid': r'CCCCCC1=CC(=C(C(=C1C(=O)O)O)CC(O2)C2(CC/C=C(C)\C)C)O',
            'gamma-eudesmyl cannabigerolate': r'CCCCCC1=CC(=C(C(=C1C(=O)OC(C)(C)C2CC3=C(C)CCCC3(C)CC2)O)C/C=C(\C)/CCC=C(C)C)O',
            'gamma-cadinyl cannabigerolate': r'CCCCCC1=CC(=C(C(=C1C(=O)OC2(C)CC(C(C)C)C3C=C(C)CCC3C2)O)C/C=C(\C)/CCC=C(C)C)O',
            'sesquicannabigerol': r'CCCCCC1=CC(=C(C(=C1)O)C/C=C(\C)/CC/C=C(\C)/CCC=C(C)C)O',
            'deprenyl O-methyl cannabigerolic acid_amorfrutin 2': 'CCCCCC1=CC(=C(C(=C1C(=O)O)O)CC=C(C)C)OC',
            '5-acetyl-4-hydroxycannabigerol': r'CCCCCC1=CC(=C(C(=C1O)OC(C)=O)C/C=C(\C)/CCC=C(C)C)O',
            'acetylcannabigeroquinol': r'CCCCCC1=CC(C(=C(C1=O)OC(C)=O)C/C=C(\C)/CCC=C(C)C)=O',
            'cannabigeroquinone': r'CCCCCC1=CC(C(=C(C1=O))C/C=C(\C)/CCC=C(C)C)=O',
            'abnormal cannabigerol': r'CCCCCC1=C(C/C=C(\C)/CCC=C(C)C)C(=CC(=C1)O)O',
            'acetyl abnormal hydrocannabigeroquinol': r'CCCCCC1=C(C/C=C(\C)/CCC=C(C)C)C(=CC(=C1OC(C)=O)O)O',
            'abnormal cannabigeroquinol': r'CCCCCC1=C(C/C=C(\C)/CCC=C(C)C)C(C=C(C1=O)O)=O',
            '2-hydroxy-1,2-dihydrocannabichromene_cyclo-CBG': r'CCCCCC1=CC(=C2CC(O)C(C)(CCC=C(C)C)OC2=C1)O',
            'cannabichromenic acid': 'CCCCCC1=C(C(O)=O)C(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
            'cannabiorcichromene': 'CC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
            'cannabiorcichromenic acid': 'CC1=C(C(O)=O)C(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
            'chlorcannabiorcichromenic acid': 'CC1=C(C(O)=O)C(=C2C=CC(OC2=C1Cl)(C)CCC=C(C)CC)O',
            '4-acetoxycannabichromene': 'CCCCCC1=CC(=C2C=CC(OC2=C1OC(C)=O)(C)CCC=C(C)C)O',
            'anthopogochromenic acid': 'O=C(O)C1=C(C=C2C=CC(OC2=C1C)(C)CCC=C(C)C)O',
            'confluentin': r'CC1=CC(=C2C=CC(OC2=C1)(C)CC/C=C(\C)/CCC=C(C)C)O',
            'daurichromenic acid': r'CC1=C(C(O)=O)C(=C2C=CC(OC2=C1)(C)CC/C=C(\C)/CCC=C(C)C)O',
            '8-hydroxyisocannabichromene': 'CCCCCC1=CC(=C2C=CC(OC2=C1)(C)CCC(O)C(=C)C)O',
            'o-methylcannabidiol': 'CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)OC',
            'o-propylcannabidiol': 'CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)OCCC',
            'o-pentylcannabidiol': 'CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)OCCCCC',
            'ferruginene C': 'CC1=CC(C(CC1)C(=C)CCC(C(=C)C)O)C2=C(C=C(C=C2O)C)O',
            'cannabioxepane': 'CCCCCC1=CC2=C3C(=C1)OCC(=C)C4=C3C(=C(C=C4)C)O2',
            '10-hydroxy delta-8-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C(O)C(=CCC3C(OC2=C1)(C)C)C)O',
            '11-acetoxy delta-8-tetrahydrocannabinolic acid': 'CCCCCC1=CC2=C(C3CC(=CCC3C(O2)(C)C)COC(C)=O)C(=C1C(=O)O)O',
            '8-hydroxy delta-9-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C=C(C(O)CC3C(OC2=C1)(C)C)C)O',
            '8-oxo delta-9-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C=C(C(=O)CC3C(OC2=C1)(C)C)C)O',
            'o-propyl delta-9-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)OCCC',
            'o-pentyl delta-9-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)OCCCCC',
            '2-formyl delta-9-tetrahydrocannabinol': 'CCCCCC1=C(C=O)C(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
            'fenchyl delta-9-tetrahydrocannabinolate': 'CCCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1C(=O)OC4C(C5CCC4(C5)C)(C)C)O',
            'bornyl delta-9-tetrahydrocannabinolate': 'CCCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1C(=O)OC4CC5CCC4(C5(C)C)C)O',
            'alpha-terpinyl delta-9-tetrahydrocannabinolate': 'CCCCCC1=C(C(OC(C)(C)C2CC=C(C)CC2)=O)C(=C3C4C=C(CCC4C(OC3=C1)(C)C)C)O',
            '4-terpinyl delta-9-tetrahydrocannabinolate': 'CCCCCC1=C(C(OC(C)C2(C)CC=C(C)CC2)=O)C(=C3C4C=C(CCC4C(OC3=C1)(C)C)C)O',
            'gamma-eudesmyl delta-9-tetrahydrocannabinolate': 'CCCCCC1=C(C(OC(C)(C)C2CC3=C(C)CCCC3(C)CC2)=O)C(=C4C5C=C(CCC5C(OC4=C1)(C)C)C)O',
            'alpha-cadinyl delta-9-tetrahydrocannabinolate': 'CCCCCC1=C(C(OC2(C)CC(C(C)C)C3C=C(C)CCC3C2)=O)C(=C4C5C=C(CCC5C(OC4=C1)(C)C)C)O',
            'tetrahydrocannabinol epoxide': 'CCCCCC1=CC(=C2C3C(O4)C4(C)CCC3C(OC2=C1)(C)(C))O',
            'hexahydrocannabinol': 'CCCCCC1=CC(=C2C3CC(CCC3C(OC2=C1)(C)C)C)O',
            'hydroxy-delta-9,11-tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3C(O)C(CCC3C(OC2=C1)(C)C)=C)O',
            'methylen-bis-delta-9-tetrahydrocannabinol_cannabisol': 'CCCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1CC4=C(C5=C(C=C4CCCCC)OC(C6C5C=C(CC6)C)(C)C)O)O',
            'cannabiorcicitran': 'CC1=CC2=C3C4CC(CCC4C(O2)(C)C)(OC3=C1)C',
            'bis-nor-cannabitriol': 'CCCC1=CC(=C2C(=C1)OC(C3=C2C(C(CC3)(O)C)O)(C)C)O',
            '10-o-ethyl-bis-nor-cannabitriol': 'CCCC1=CC(=C2C(=C1)OC(C3=C2C(C(CC3)(O)C)OCC)(C)C)O',
            'isocannabitriol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C(C(C(O)C3)(O)C))(C)C)O',
            '10-o-ethyl-cannabitriol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C(C(CC3)(O)C)OCC)(C)C)O',
            '9,10-anhydrocannabitriol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C4O(C4(CC3)C))(C)C)O',
            '7,8-dehydro-10-o-ethylcannabitriol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C(C(C=C3)(O)C)OCC)(C)C)O',
            'delta-7-isotetrahydrocannabivarin': 'CCCC1=CC(=C2C3CC(CCC3C(=C)C)(OC2=C1)C)O',
            'delta-7-isotetrahydrocannabinol': 'CCCCCC1=CC(=C2C3CC(CCC3C(=C)C)(OC2=C1)C)O',
            'cannabiglendol': 'CCCC1=CC(=C2C3CC(CCC3C(C)(C)O)(OC2=C1)C)O',
            'bis-nor-cannabielsoin': 'CCCC1=CC(=C2C3C(CCC(C3OC2=C1)(C)O)C(=C)C)O',
            'bis-nor-cannabielsoic acid B': 'CCCC1=CC(=C2C3C(CCC(C3OC2=C1C(=O)O)(C)O)C(=C)C)O',
            'ferruginene A': 'CC1=CC(=C2C3C(CCC(C3OC2=C1)(C)O)C(=C)CCC(C(=C)C)O)O',
            'ferruginene B': 'CC1=CC(=C2C3C(CCC(C3OC2=C1)(C)O)C(=C)CC=CC(C)(C)O)O',
            'cannabiorcicyclol': 'CC1=CC(=C2C3C4C(C3(C)C)CCC4(OC2=C1)C)O',
            'cannabiorcicyclolic acid': 'CC1=CC2=C(C3C4C(C3(C)C)CCC4(O2)C)C(=C1C(=O)O)O',
            'anthopogocyclolic acid': 'OC1=CC(=C2C3C4C(C3(C)C)CCC4(OC2=C1)C)C',
            'rhododaurichromanic acid A': 'OC1=C(C(O)=O)C(=C2C3C4C(C3(CC=C(C)C)C)CCC4(OC2=C1)C)C',
            'nor-cannabivarin': 'CCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)O',
            'o-methylcannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)OC',
            'o-propylcannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)OCCC',
            'o-pentylcannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)OCCCCC',
            '7-hydroxycannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3(O))C)(C)C)O',
            '8-hydroxycannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C(O)=C3)C)(C)C)O',
            '8-hydroxycannabinolic acid': 'CCCCCC1=C(C(O)=O)C(=C2C(=C1)OC(C3=C2C=C(C(O)=C3)C)(C)C)O',
            '7,8-dihydrocannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(CC3)C)(C)C)O',
            '4-terpenylcannabinolate': 'CCCCCC1=CC2=C(C3=C(C=CC(=C3)C)C(O2)(C)C)C(=C1C(=O)OC4(CCC(=CC4)C)C(C)C)O',
            'cannabichromanone': 'CCCCCC1=CC(=C2C(=C1)OC(C(C2=O)CCC(=O)C)(C)C)O',
            'bis-nor-cannabichromanone': 'CCCC1=CC(=C2C(=C1)OC(C(C2=O)CCC(=O)C)(C)C)O',
            'cannabichromanone B': 'CCCCCC1=CC(=C2C(=C1)OC(C(C2=O)C(O)CC(=O)C)(C)C)O',
            'cannabichromanone C': 'CCCCCC1=CC(=C2C(=C1)OC(C(C2=O)CC(=O)C(=O)C)(C)C)O',
            'cannabicoumaronone': 'CCCCCC1=CC2=C3C(=C1)OC(C(C3=CO2)CCC(=O)C)(C)C',
            'cannabicoumarononic acid': 'CCCCCC1=C(C(O)=O)C2=C3C(=C1)OC(C(C3=CO2)CCC(=O)C)(C)C',
            'cannabimovone': r'CCCCCC1=CC(=C(C(=C1)O)C2C(CC(C2O)C(=O)C)C(=C)C)O',
            'anhydrocannabimovone': 'CCCCCC1=CC(=C2C(=C1)OC3C2C(CC3C(=O)C)C(=C)C)O',
            'amorfrutin 1_amorfrutin A': 'CC(=CCC1=C(C=C(C(=C1O)C(=O)O)CCC2=CC=CC=C2)OC)C',
            'demethyldecarboxyamorfrutin A': 'CC(=CCC1=C(C=C(C=C1O)CCC2=CC=CC=C2)O)C',
            'demethylamorfrutin A': 'CC(=CCC1=C(C(C(O)=O)=C(C=C1O)CCC2=CC=CC=C2)O)C',
            'decarboxyamorfrutin A': 'CC(=CCC1=C(C=C(C=C1OC)CCC2=CC=CC=C2)O)C',
            'amorfrutin 4_amorfrutin C': 'CC(=CCC1=C(C(=C(C(=C1OC)CC=C(C)C)O)C(=O)OC)CCC2=CC=CC=C2)C',
            'heli-cannabigerol': r'CC(=CCC/C(=C\CC1=C(C=C(C=C1O)CCC2=CC=CC=C2)O)C)C',
            'hydroxy-heli-cannabigerol': r'CC(=CCC/C(=C\CC1=C(C=C(C=C1O)CCC2=CC=C(O)C=C2)O)C)C',
            'amorfrutin B': r'CC(=CCC/C(=C/CC1=C(C=C(C(=C1O)C(=O)O)CCC2=CC=CC=C2)OC)/C)C',
            'demethylamorfrutin B': r'CC(=CCC/C(=C\CC1=C(C(C(O)=O)=C(C=C1O)CCC2=CC=CC=C2)O)C)C',
            'decarboxyamorfrutin B': r'CC(=CCC/C(=C\CC1=C(C=C(C=C1OC)CCC2=CC=CC=C2)O)C)C',
            'chiricanin A': 'CC(=CCC1=C(C=C(C=C1O)C=CC2=CC=CC=C2)O)C',
            'arachidin-2': 'CC(=CCC1=C(C=C(C=C1O)C=CC2=CC=C(C=C2)O)O)C',
            'glepidotin C': r'CC(=C)C(CC1=C(C=C(C=C1O)CCC2=CC=CC=C2)O)O',
            'amorfrutin 3': 'CC(=C)C(CC1=C(C=C(C(=C1O)C(=O)O)CCC2=CC=CC=C2)OC)O',
            'arachidin-3': r'CC(C)/C=C/C1=C(C=C(C=C1O)/C=C/C2=CC=C(C=C2)O)O',
            'arachidin-1': r'CC(C)/C=C/C1=C(C=C(C=C1O)/C=C/C2=CC(=C(C=C2)O)O)O',
            'arachidin-4': r'CC(C)(CCC1=C(C=C(C=C1O)/C=C/C2=CC=C(C=C2)O)O)O',
            'amorfrutin D': r'CC(=C(CC/C(=C/CC1=C(C=C(C(=C1O)C(=O)O)CCC2=CC=CC=C2)OC)/C)O)C',
            'machaeridol A': r'CC1CCC(C(C1)C2=C(C=C(C=C2O)/C=C/C3=CC=CC=C3)O)C(=C)C',
            'machaeridol B': r'CC1CCC(C(C1)C2=C(C=C(C=C2O)/C=C/C3=CC=CC=C3O)O)C(=C)C',
            'machaeridol C': 'CC1CCC(C(C1)C2=C(C=C(C=C2O)C3=CC4=CC=CC=C4O3)O)C(=C)C',
            'machaeriol A': r'CC1CCC2C(C1)C3=C(C=C(C=C3OC2(C)C)/C=C/C4=CC=CC=C4)O',
            'machaeriol B': 'CC1CCC2C(C1)C3=C(C=C(C=C3OC2(C)C)C4=CC5=CC=CC=C5O4)O',
            'machaeriol C': r'CC1CCC2C(C1)C3=C(C=C(C=C3OC2(C)C)/C=C/C4=CC=CC=C4O)O',
            'machaeriol D': 'CC1CC2C(CC1O)C(OC3=CC(=CC(=C23)O)C4=CC5=CC=CC=C5O4)(C)C',
            'tetrahydrocannabiphorol': 'CCCCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
            'cannabidiphorol': 'CCCCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O'
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'cannabinerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]/[#6]=[#6](/[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'cannabinerolic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6]/[#6]=[#6](/[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'carmagerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6](-[#6](-[#6])(-[#6])-[#8])-[#8])-[#8]',
            'rac-6-epoxycannabigerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]-[#6]1-[#8]-[#6]-1(-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#6])-[#8]',
            'rac-6-epoxycannabigerolic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6]-[#6]1-[#8]-[#6]-1(-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#6])-[#8]',
            'rac-6-epoxycannabinerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]-[#6]1-[#8]-[#6]-1(-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#6])-[#8]',
            'rac-6-epoxycannabinerolic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6]-[#6]1-[#8]-[#6]-1(-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#6])-[#8]',
            'gamma-eudesmyl cannabigerolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])-[#6]1-[#6]-[#6]2=[#6](-[#6])-[#6]-[#6]-[#6]-[#6]-2(-[#6])-[#6]-[#6]-1)-[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'gamma-cadinyl cannabigerolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#6](=[#8])-[#8]-[#6]1(-[#6])-[#6]-[#6](-[#6](-[#6])-[#6])-[#6]2-[#6]=[#6](-[#6])-[#6]-[#6]-[#6]-2-[#6]-1)-[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'sesquicannabigerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'deprenyl o-methyl cannabigerolic acid_amorfrutin 2': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]-[#6]',
            '5-acetyl-4-hydroxycannabigerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1-[#8])-[#8]-[#6](-[#6])=[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'acetylcannabigeroquinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1=[#6]-[#6](-[#6](=[#6](-[#6]-1=[#8])-[#8]-[#6](-[#6])=[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])=[#8]',
            'cannabigeroquinone': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1=[#6]-[#6](-[#6](=[#6]-[#6]-1=[#8])-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])=[#8]',
            'abnormal cannabigerol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6]):[#6](:[#6]:[#6](:[#6]:1)-[#8])-[#8]',
            'acetyl abnormal hydrocannabigeroquinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6]):[#6](:[#6]:[#6](:[#6]:1-[#8]-[#6](-[#6])=[#8])-[#8])-[#8]',
            'abnormal cannabigeroquinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1=[#6](-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#6](-[#6]=[#6](-[#6]-1=[#8])-[#8])=[#8]',
            '2-hydroxy-1,2-dihydrocannabichromene_cyclo-cbg': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]-[#6](-[#8])-[#6](-[#6])(-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]-[#6]:2:[#6]:1)-[#8]',
            'cannabichromenic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'cannabiorcichromene': '[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'cannabiorcichromenic acid': '[#6]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'chlorcannabiorcichromenic acid': '[#6]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1-[#17])(-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6]-[#6])-[#8]',
            '4-acetoxycannabichromene': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1-[#8]-[#6](-[#6])=[#8])(-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'anthopogochromenic acid': '[#8]=[#6](-[#8])-[#6]1:[#6](:[#6]:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1-[#6])(-[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'confluentin': '[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6]-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            'daurichromenic acid': '[#6]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6]-[#6]/[#6]=[#6](\[#6])-[#6]-[#6]-[#6]=[#6](-[#6])-[#6])-[#8]',
            '8-hydroxyisocannabichromene': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]=[#6]-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6]-[#6]-[#6](-[#8])-[#6](=[#6])-[#6])-[#8]',
            'o-methylcannabidiol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]1-[#6]=[#6](-[#6]-[#6]-[#6]-1-[#6](=[#6])-[#6])-[#6])-[#8]-[#6]',
            'o-propylcannabidiol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]1-[#6]=[#6](-[#6]-[#6]-[#6]-1-[#6](=[#6])-[#6])-[#6])-[#8]-[#6]-[#6]-[#6]',
            'o-pentylcannabidiol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]1-[#6]=[#6](-[#6]-[#6]-[#6]-1-[#6](=[#6])-[#6])-[#6])-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]',
            'ferruginene c': '[#6]-[#6]1=[#6]-[#6](-[#6](-[#6]-[#6]-1)-[#6](=[#6])-[#6]-[#6]-[#6](-[#6](=[#6])-[#6])-[#8])-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6])-[#8]',
            'cannabioxepane': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6]3:[#6](:[#6]:1)-[#8]-[#6]-[#6](=[#6])-[#6]1:[#6]:3:[#6](:[#6](:[#6]:[#6]:1)-[#6]):[#8]:2',
            '10-hydroxy delta-8-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6](-[#8])-[#6](=[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            '11-acetoxy delta-8-tetrahydrocannabinolic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6](-[#6]3-[#6]-[#6](=[#6]-[#6]-[#6]-3-[#6](-[#8]-2)(-[#6])-[#6])-[#6]-[#8]-[#6](-[#6])=[#8]):[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8]',
            '8-hydroxy delta-9-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6](-[#8])-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            '8-oxo delta-9-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6](=[#8])-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'o-propyl delta-9-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]-[#6]-[#6]-[#6]',
            'o-pentyl delta-9-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '2-formyl delta-9-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6]=[#8]):[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'fenchyl delta-9-tetrahydrocannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6](-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-2)(-[#6])-[#6])-[#6]):[#6](:[#6]:1-[#6](=[#8])-[#8]-[#6]1-[#6](-[#6]2-[#6]-[#6]-[#6]-1(-[#6]-2)-[#6])(-[#6])-[#6])-[#8]',
            'bornyl delta-9-tetrahydrocannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6](-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-2)(-[#6])-[#6])-[#6]):[#6](:[#6]:1-[#6](=[#8])-[#8]-[#6]1-[#6]-[#6]2-[#6]-[#6]-[#6]-1(-[#6]-2(-[#6])-[#6])-[#6])-[#8]',
            'alpha-terpinyl delta-9-tetrahydrocannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6]2-[#6]-[#6]=[#6](-[#6])-[#6]-[#6]-2)=[#8]):[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            '4-terpinyl delta-9-tetrahydrocannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8]-[#6](-[#6])-[#6]2(-[#6])-[#6]-[#6]=[#6](-[#6])-[#6]-[#6]-2)=[#8]):[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'gamma-eudesmyl delta-9-tetrahydrocannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6]2-[#6]-[#6]3=[#6](-[#6])-[#6]-[#6]-[#6]-[#6]-3(-[#6])-[#6]-[#6]-2)=[#8]):[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'alpha-cadinyl delta-9-tetrahydrocannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8]-[#6]2(-[#6])-[#6]-[#6](-[#6](-[#6])-[#6])-[#6]3-[#6]=[#6](-[#6])-[#6]-[#6]-[#6]-3-[#6]-2)=[#8]):[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'tetrahydrocannabinol epoxide': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]4-[#8]-[#6]-4(-[#6])-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#8]',
            'hexahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]-[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'hydroxy-delta-9,11-tetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6](-[#8])-[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])=[#6])-[#8]',
            'methylen-bis-delta-9-tetrahydrocannabinol_cannabisol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6](-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-2)(-[#6])-[#6])-[#6]):[#6](:[#6]:1-[#6]-[#6]1:[#6](:[#6]2:[#6](:[#6]:[#6]:1-[#6]-[#6]-[#6]-[#6]-[#6])-[#8]-[#6](-[#6]1-[#6]-2-[#6]=[#6](-[#6]-[#6]-1)-[#6])(-[#6])-[#6])-[#8])-[#8]',
            'cannabiorcicitran': '[#6]-[#6]1:[#6]:[#6]2:[#6]3-[#6]4-[#6]-[#6](-[#6]-[#6]-[#6]-4-[#6](-[#8]-2)(-[#6])-[#6])(-[#8]-[#6]:3:[#6]:1)-[#6]',
            'bis-nor-cannabitriol': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6](-[#6](-[#6]-[#6]-1)(-[#8])-[#6])-[#8])(-[#6])-[#6])-[#8]',
            '10-o-ethyl-bis-nor-cannabitriol': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6](-[#6](-[#6]-[#6]-1)(-[#8])-[#6])-[#8]-[#6]-[#6])(-[#6])-[#6])-[#8]',
            'isocannabitriol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6]-[#6](-[#6](-[#8])-[#6]-1)(-[#8])-[#6])(-[#6])-[#6])-[#8]',
            '10-o-ethyl-cannabitriol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6](-[#6](-[#6]-[#6]-1)(-[#8])-[#6])-[#8]-[#6]-[#6])(-[#6])-[#6])-[#8]',
            '9,10-anhydrocannabitriol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6]2-[#8]-[#6]-2(-[#6]-[#6]-1)-[#6])(-[#6])-[#6])-[#8]',
            '7,8-dehydro-10-o-ethylcannabitriol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6](-[#6](-[#6]=[#6]-1)(-[#8])-[#6])-[#8]-[#6]-[#6])(-[#6])-[#6])-[#8]',
            'delta-7-isotetrahydrocannabivarin': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]-[#6](-[#6]-[#6]-[#6]-3-[#6](=[#6])-[#6])(-[#8]-[#6]:2:[#6]:1)-[#6])-[#8]',
            'delta-7-isotetrahydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]-[#6](-[#6]-[#6]-[#6]-3-[#6](=[#6])-[#6])(-[#8]-[#6]:2:[#6]:1)-[#6])-[#8]',
            'cannabiglendol': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]-[#6](-[#6]-[#6]-[#6]-3-[#6](-[#6])(-[#6])-[#8])(-[#8]-[#6]:2:[#6]:1)-[#6])-[#8]',
            'bis-nor-cannabielsoin': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6](-[#6]-[#6]-[#6](-[#6]-3-[#8]-[#6]:2:[#6]:1)(-[#6])-[#8])-[#6](=[#6])-[#6])-[#8]',
            'bis-nor-cannabielsoic acid b': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6](-[#6]-[#6]-[#6](-[#6]-3-[#8]-[#6]:2:[#6]:1-[#6](=[#8])-[#8])(-[#6])-[#8])-[#6](=[#6])-[#6])-[#8]',
            'ferruginene a': '[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6](-[#6]-[#6]-[#6](-[#6]-3-[#8]-[#6]:2:[#6]:1)(-[#6])-[#8])-[#6](=[#6])-[#6]-[#6]-[#6](-[#6](=[#6])-[#6])-[#8])-[#8]',
            'ferruginene b': '[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6](-[#6]-[#6]-[#6](-[#6]-3-[#8]-[#6]:2:[#6]:1)(-[#6])-[#8])-[#6](=[#6])-[#6]-[#6]=[#6]-[#6](-[#6])(-[#6])-[#8])-[#8]',
            'cannabiorcicyclol': '[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]4-[#6](-[#6]-3(-[#6])-[#6])-[#6]-[#6]-[#6]-4(-[#8]-[#6]:2:[#6]:1)-[#6])-[#8]',
            'cannabiorcicyclolic acid': '[#6]-[#6]1:[#6]:[#6]2:[#6](-[#6]3-[#6]4-[#6](-[#6]-3(-[#6])-[#6])-[#6]-[#6]-[#6]-4(-[#8]-2)-[#6]):[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8]',
            'anthopogocyclolic acid': '[#8]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]4-[#6](-[#6]-3(-[#6])-[#6])-[#6]-[#6]-[#6]-4(-[#8]-[#6]:2:[#6]:1)-[#6])-[#6]',
            'rhododaurichromanic acid a': '[#8]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]2-[#6]3-[#6]4-[#6](-[#6]-3(-[#6]-[#6]=[#6](-[#6])-[#6])-[#6])-[#6]-[#6]-[#6]-4(-[#8]-[#6]:2:[#6]:1)-[#6])-[#6]',
            'nor-cannabivarin': '[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6]:[#6]:1)-[#6])(-[#6])-[#6])-[#8]',
            'o-methylcannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6]:[#6]:1)-[#6])(-[#6])-[#6])-[#8]-[#6]',
            'o-propylcannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6]:[#6]:1)-[#6])(-[#6])-[#6])-[#8]-[#6]-[#6]-[#6]',
            'o-pentylcannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6]:[#6]:1)-[#6])(-[#6])-[#6])-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]',
            '7-hydroxycannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6])(-[#6])-[#6])-[#8]',
            '8-hydroxycannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6](-[#8]):[#6]:1)-[#6])(-[#6])-[#6])-[#8]',
            '8-hydroxycannabinolic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1:[#6]-2:[#6]:[#6](:[#6](-[#8]):[#6]:1)-[#6])(-[#6])-[#6])-[#8]',
            '7,8-dihydrocannabinol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6]1=[#6]-2-[#6]=[#6](-[#6]-[#6]-1)-[#6])(-[#6])-[#6])-[#8]',
            '4-terpenylcannabinolate': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6](-[#6]3:[#6](:[#6]:[#6]:[#6](:[#6]:3)-[#6])-[#6](-[#8]-2)(-[#6])-[#6]):[#6](:[#6]:1-[#6](=[#8])-[#8]-[#6]1(-[#6]-[#6]-[#6](=[#6]-[#6]-1)-[#6])-[#6](-[#6])-[#6])-[#8]',
            'cannabichromanone': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6](-[#6]-2=[#8])-[#6]-[#6]-[#6](=[#8])-[#6])(-[#6])-[#6])-[#8]',
            'bis-nor-cannabichromanone': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6](-[#6]-2=[#8])-[#6]-[#6]-[#6](=[#8])-[#6])(-[#6])-[#6])-[#8]',
            'cannabichromanone b': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6](-[#6]-2=[#8])-[#6](-[#8])-[#6]-[#6](=[#8])-[#6])(-[#6])-[#6])-[#8]',
            'cannabichromanone c': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6](-[#6](-[#6]-2=[#8])-[#6]-[#6](=[#8])-[#6](=[#8])-[#6])(-[#6])-[#6])-[#8]',
            'cannabicoumaronone': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]2:[#6]3:[#6](:[#6]:1)-[#8]-[#6](-[#6](-[#6]:3:[#6]:[#8]:2)-[#6]-[#6]-[#6](=[#8])-[#6])(-[#6])-[#6]',
            'cannabicoumarononic acid': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6](-[#6](-[#8])=[#8]):[#6]2:[#6]3:[#6](:[#6]:1)-[#8]-[#6](-[#6](-[#6]:3:[#6]:[#8]:2)-[#6]-[#6]-[#6](=[#8])-[#6])(-[#6])-[#6]',
            'cannabimovone': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]1-[#6](-[#6]-[#6](-[#6]-1-[#8])-[#6](=[#8])-[#6])-[#6](=[#6])-[#6])-[#8]',
            'anhydrocannabimovone': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2:[#6](:[#6]:1)-[#8]-[#6]1-[#6]-2-[#6](-[#6]-[#6]-1-[#6](=[#8])-[#6])-[#6](=[#6])-[#6])-[#8]',
            'amorfrutin 1_amorfrutin a': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6](:[#6]:1-[#8])-[#6](=[#8])-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8]-[#6])-[#6]',
            'demethyldecarboxyamorfrutin a': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6]',
            'demethylamorfrutin a': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]:[#6]:1-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6]',
            'decarboxyamorfrutin a': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8]-[#6])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6]',
            'amorfrutin 4_amorfrutin c': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6](:[#6](:[#6](:[#6]:1-[#8]-[#6])-[#6]-[#6]=[#6](-[#6])-[#6])-[#8])-[#6](=[#8])-[#8]-[#6])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]',
            'heli-cannabigerol': '[#6]-[#6](=[#6]-[#6]-[#6]/[#6](=[#6]\[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6])-[#6]',
            'hydroxy-heli-cannabigerol': '[#6]-[#6](=[#6]-[#6]-[#6]/[#6](=[#6]\[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]):[#6]:[#6]:1)-[#8])-[#6])-[#6]',
            'amorfrutin b': '[#6]-[#6](=[#6]-[#6]-[#6]/[#6](=[#6]/[#6]-[#6]1:[#6](:[#6]:[#6](:[#6](:[#6]:1-[#8])-[#6](=[#8])-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8]-[#6])-[#6])-[#6]',
            'demethylamorfrutin b': '[#6]-[#6](=[#6]-[#6]-[#6]/[#6](=[#6]\[#6]-[#6]1:[#6](:[#6](-[#6](-[#8])=[#8]):[#6](:[#6]:[#6]:1-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6])-[#6]',
            'decarboxyamorfrutin b': '[#6]-[#6](=[#6]-[#6]-[#6]/[#6](=[#6]\[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8]-[#6])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6])-[#6]',
            'chiricanin a': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]=[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6]',
            'arachidin-2': '[#6]-[#6](=[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]=[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8])-[#8])-[#6]',
            'glepidotin c': '[#6]-[#6](=[#6])-[#6](-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#8]',
            'amorfrutin 3': '[#6]-[#6](=[#6])-[#6](-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6](:[#6]:1-[#8])-[#6](=[#8])-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8]-[#6])-[#8]',
            'arachidin-3': '[#6]-[#6](-[#6])/[#6]=[#6]/[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])/[#6]=[#6]/[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8])-[#8]',
            'arachidin-1': '[#6]-[#6](-[#6])/[#6]=[#6]/[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])/[#6]=[#6]/[#6]1:[#6]:[#6](:[#6](:[#6]:[#6]:1)-[#8])-[#8])-[#8]',
            'arachidin-4': '[#6]-[#6](-[#6])(-[#6]-[#6]-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])/[#6]=[#6]/[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8])-[#8])-[#8]',
            'amorfrutin d': '[#6]-[#6](=[#6](-[#6]-[#6]/[#6](=[#6]/[#6]-[#6]1:[#6](:[#6]:[#6](:[#6](:[#6]:1-[#8])-[#6](=[#8])-[#8])-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8]-[#6])-[#6])-[#8])-[#6]',
            'machaeridol a': '[#6]-[#6]1-[#6]-[#6]-[#6](-[#6](-[#6]-1)-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])/[#6]=[#6]/[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8])-[#6](=[#6])-[#6]',
            'machaeridol b': '[#6]-[#6]1-[#6]-[#6]-[#6](-[#6](-[#6]-1)-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])/[#6]=[#6]/[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8])-[#8])-[#6](=[#6])-[#6]',
            'machaeridol c': '[#6]-[#6]1-[#6]-[#6]-[#6](-[#6](-[#6]-1)-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8])-[#6]1:[#6]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#8]:1)-[#8])-[#6](=[#6])-[#6]',
            'machaeriol a': '[#6]-[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]-1)-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8]-[#6]-2(-[#6])-[#6])/[#6]=[#6]/[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#8]',
            'machaeriol b': '[#6]-[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]-1)-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8]-[#6]-2(-[#6])-[#6])-[#6]1:[#6]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#8]:1)-[#8]',
            'machaeriol c': '[#6]-[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]-1)-[#6]1:[#6](:[#6]:[#6](:[#6]:[#6]:1-[#8]-[#6]-2(-[#6])-[#6])/[#6]=[#6]/[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8])-[#8]',
            'machaeriol d': '[#6]-[#6]1-[#6]-[#6]2-[#6](-[#6]-[#6]-1-[#8])-[#6](-[#8]-[#6]1:[#6]:[#6](:[#6]:[#6](:[#6]:1-2)-[#8])-[#6]1:[#6]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#8]:1)(-[#6])-[#6]',
            'tetrahydrocannabiphorol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6]2-[#6]3-[#6]=[#6](-[#6]-[#6]-[#6]-3-[#6](-[#8]-[#6]:2:[#6]:1)(-[#6])-[#6])-[#6])-[#8]',
            'cannabidiphorol': '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6](:[#6](:[#6](:[#6]:1)-[#8])-[#6]1-[#6]=[#6](-[#6]-[#6]-[#6]-1-[#6](=[#6])-[#6])-[#6])-[#8]',
        }

        return smarts