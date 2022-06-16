#!/usr/bin/env python3
#
# GlobalChem - FDA List Five
#
# -------------------------

class FDAListFive(object):

    def __init__(self):

        self.name = 'fda_list_five'

    @staticmethod
    def get_smiles():

        smiles = {
            'alumina': 'O=[Al]O[Al]=O',
            'aluminum powder': '[Al]',
            'annatto extract': 'CC(=CC=CC=C(C)C=CC=C(C)C=CC(=O)O)C=CC=C(C)C=CC(=O)O',
            'beta-carotene': 'CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2=C(CCCC2(C)C)C)C)C',
            'bismuth oxychloride': 'O=[Bi].Cl',
            'bronze powder': '[Cu].[Sn]',
            'calcium carbonate': 'C(=O)([O-])[O-].[Ca+2]',
            'canthaxanthin': 'CC1=C(C(CCC1=O)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2=C(C(=O)CCC2(C)C)C)C)C',
            'caramel': '*',
            'carmine': 'CC1=C2C(=CC(=C1C(=O)O)O)C(=O)C3=C(C2=O)C(=C(C(=C3O)O)C4C(C(C(C(O4)CO)O)O)O)O',
            'chlorophyllin, copper complex': 'CCC1=C(C2=NC1=CC3=C(C(=C([O-])[O-])C(=N3)C(=C4C(C(C(=CC5=NC(=C2)C(=C5C)C=C)N4)C)CCC(=O)O)CC(=O)O)C)C.[Cu+2]',
            'chromium hydroxide, green': 'O.O.[O-2].[O-2].[O-2].[Cr+3].[Cr+3]',
            'chromium oxides greens': 'O=[Cr]O[Cr]=O',
            'cochineal extract': 'CC1=C2C(=CC(=C1C(=O)O)O)C(=O)C3=C(C2=O)C(=C(C(=C3O)O)C4C(C(C(C(O4)CO)O)O)O)O',
            'copper, metallic powder': '[Cu]',
            'potassium sodium copper chlorophyllin': 'CCC1=C(C2=NC1=CC3=C(C(=C([N-]3)C(=C4C(C(C(=N4)C=C5C(=C(C(=C2)[N-]5)C=C)C)C)CCC(=O)[O-])CC(=O)[O-])C(=O)[O-])C)C.[Na+].[Na+].[Na+].[Cu+2]',
            'dihydroxyacetone': 'C(C(=O)CO)O',
            'ferric ammonium ferrocyanide': '[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[NH4+].[Fe+2].[Fe+3]',
            'ferric ferrocyanide': '[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[Fe+3].[Fe+3]',
            'guanine': 'C1=NC2=C(N1)C(=O)NC(=N2)N',
            'mica': '[O-2].O=[Al]O[Al]=O.O=[Si]=O.[K+].[K+]',
            'mica-based pearlescent pigment': '*',
            'pyrophyllite': 'O[Si](=O)[O-].O[Si](=O)[O-].[Al+2]',
            'synthetic iron oxide': 'O=[Fe]O[Fe]=O',
            'talc': 'O.O=[Mg].O=[Mg].O=[Mg].O=[Si]=O.O=[Si]=O.O=[Si]=O.O=[Si]=O',
            'titanium dioxide': 'O=[Ti]=O',
            'zinc oxide': 'O=[Zn]',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'alumina': '[#8]=[Al]-[#8]-[Al]=[#8]',
            'aluminum powder': '[Al]',
            'annatto extract': '[#6]-[#6](=[#6]-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6](=[#8])-[#8])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6](=[#8])-[#8]',
            'beta-carotene': '[#6]-[#6]1=[#6](-[#6](-[#6]-[#6]-[#6]-1)(-[#6])-[#6])-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]1=[#6](-[#6]-[#6]-[#6]-[#6]-1(-[#6])-[#6])-[#6])-[#6])-[#6]',
            'bismuth oxychloride': '[#8]=[Bi].[#17]',
            'bronze powder': '[Cu].[Sn]',
            'calcium carbonate': '[#6](=[#8])(-[#8-])-[#8-].[Ca+2]',
            'canthaxanthin': '[#6]-[#6]1=[#6](-[#6](-[#6]-[#6]-[#6]-1=[#8])(-[#6])-[#6])-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]1=[#6](-[#6](=[#8])-[#6]-[#6]-[#6]-1(-[#6])-[#6])-[#6])-[#6])-[#6]',
            'caramel': '[#0]',
            'carmine': '[#6]-[#6]1:[#6]2:[#6](:[#6]:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6](=[#8])-[#6]1:[#6](-[#6]-2=[#8]):[#6](:[#6](:[#6](:[#6]:1-[#8])-[#8])-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8])-[#8]',
            'chlorophyllin, copper complex': '[#6]-[#6]-[#6]1=[#6](-[#6]2=[#7]-[#6]-1=[#6]-[#6]1=[#6](-[#6](=[#6](-[#8-])-[#8-])-[#6](=[#7]-1)-[#6](=[#6]1-[#6](-[#6](-[#6](=[#6]-[#6]3=[#7]-[#6](=[#6]-2)-[#6](=[#6]-3-[#6])-[#6]=[#6])-[#7]-1)-[#6])-[#6]-[#6]-[#6](=[#8])-[#8])-[#6]-[#6](=[#8])-[#8])-[#6])-[#6].[Cu+2]',
            'chromium hydroxide, green': '[#8].[#8].[#8-2].[#8-2].[#8-2].[Cr+3].[Cr+3]',
            'chromium oxides greens': '[#8]=[Cr]-[#8]-[Cr]=[#8]',
            'cochineal extract': '[#6]-[#6]1:[#6]2:[#6](:[#6]:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6](=[#8])-[#6]1:[#6](-[#6]-2=[#8]):[#6](:[#6](:[#6](:[#6]:1-[#8])-[#8])-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8])-[#8]',
            'copper, metallic powder': '[Cu]',
            'potassium sodium copper chlorophyllin': '[#6]-[#6]-[#6]1=[#6](-[#6]2:[#7]:[#6]-1:[#6]:[#6]1:[#6](:[#6](:[#6](:[#7-]:1):[#6](:[#6]1-[#6](-[#6](-[#6](:[#7]:1):[#6]:[#6]1:[#6](:[#6](:[#6](:[#6]:2):[#7-]:1)-[#6]=[#6])-[#6])-[#6])-[#6]-[#6]-[#6](=[#8])-[#8-])-[#6]-[#6](=[#8])-[#8-])-[#6](=[#8])-[#8-])-[#6])-[#6].[Na+].[Na+].[Na+].[Cu+2]',
            'dihydroxyacetone': '[#6](-[#6](=[#8])-[#6]-[#8])-[#8]',
            'ferric ammonium ferrocyanide': '[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#7H4+].[Fe+2].[Fe+3]',
            'ferric ferrocyanide': '[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[Fe+3].[Fe+3]',
            'guanine': '[#6]1:[#7]:[#6]2:[#6](:[#7H]:1):[#6](=[#8]):[#7H]:[#6](:[#7]:2)-[#7]',
            'mica': '[#8-2].[#8]=[Al]-[#8]-[Al]=[#8].[#8]=[Si]=[#8].[K+].[K+]',
            'mica-based pearlescent pigment': '[#0]',
            'pyrophyllite': '[#8]-[Si](=[#8])-[#8-].[#8]-[Si](=[#8])-[#8-].[Al+2]',
            'synthetic iron oxide': '[#8]=[Fe]-[#8]-[Fe]=[#8]',
            'talc': '[#8].[#8]=[Mg].[#8]=[Mg].[#8]=[Mg].[#8]=[Si]=[#8].[#8]=[Si]=[#8].[#8]=[Si]=[#8].[#8]=[Si]=[#8]',
            'titanium dioxide': '[#8]=[Ti]=[#8]',
            'zinc oxide': '[#8]=[Zn]',
        }

        return smarts