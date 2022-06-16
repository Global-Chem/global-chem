#!/usr/bin/env python3
#
# GlobalChem - FDA List Six
#
# -------------------------

class FDAListSix(object):

    def __init__(self):

        self.name = 'fda_list_six'

    @staticmethod
    def get_smiles():

        smiles = {
            'aluminum powder': '[Al]',
            'annatto': 'CC(=CC=CC=C(C)C=CC=C(C)C=CC(=O)O)C=CC=C(C)C=CC(=O)O',
            'bismuth citrate': '',
            'bismuth oxychloride': 'O=[Bi].Cl',
            'Bronze powder': '[Cu].[Sn]',
            'caramel': '*',
            'carmine': 'CC1=C2C(=CC(=C1C(=O)O)O)C(=O)C3=C(C2=O)C(=C(C(=C3O)O)C4C(C(C(C(O4)CO)O)O)O)O',
            'beta-carotene': 'CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2=C(CCCC2(C)C)C)C)C',
            'chromium hydroxide green': 'O.O.[O-2].[O-2].[O-2].[Cr+3].[Cr+3]',
            'chromium oxide greens': '[O-2].[O-2].[O-2].[Cr+3].[Cr+3]',
            'copper, metallic powder ': '[Cu]',
            'dihydroxyacetone': 'C(C(=O)CO)O',
            'disodium edta-copper': 'C(CN(CC(=O)[O-])CC(=O)[O-])N(CC(=O)[O-])CC(=O)[O-].[Na+].[Na+].[Cu+2]',
            'ferric ammonium ferrocyanide': '[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[NH4+].[Fe+2].[Fe+3]',
            'ferric ferrocyanide': '[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[Fe+3].[Fe+3]',
            'guaiazulene': 'CC1=C2C=CC(=C2C=C(C=C1)C(C)C)C',
            'guanine': 'C1=NC2=C(N1)C(=O)NC(=N2)N',
            'henna': 'C1=CC=C2C(=C1)C(=CC(=O)C2=O)O',
            'iron oxides': 'O=[Fe]O[Fe]=O',
            'luminescent zinc sulfide': 'S=[Zn]',
            'manganese violet': '[NH4+].[O-]P(=O)([O-])OP(=O)([O-])[O-].[Mn+3]',
            'mica': '[O-2].O=[Al]O[Al]=O.O=[Si]=O.[K+].[K+]',
            'potassium sodium copper chlorophyllin': 'CCC1=C(C2=NC1=CC3=C(C(=C([N-]3)C(=C4C(C(C(=N4)C=C5C(=C(C(=C2)[N-]5)C=C)C)C)CCC(=O)[O-])CC(=O)[O-])C(=O)[O-])C)C.[Na+].[Na+].[Na+].[Cu+2]',
            'pyrophyllite': 'O[Si](=O)[O-].O[Si](=O)[O-].[Al+2]',
            'silver': '[Ag]',
            'silver nitrate': '[N+](=O)([O-])[O-].[Ag+]',
            'titanium dioxide': 'O=[Ti]=O',
            'ultramarines': '[O-][Si]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Al+3].[Al+3].[Al+3].[Al+3].[Al+3].[Al+3].[S-]S[S-]',
            'zinc oxide': 'O=[Zn]',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'aluminum powder': '[Al]',
            'annatto': '[#6]-[#6](=[#6]-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6](=[#8])-[#8])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6](=[#8])-[#8]',
            'bismuth citrate': '',
            'bismuth oxychloride': '[#8]=[Bi].[#17]',
            'bronze powder': '[Cu].[Sn]',
            'caramel': '[#0]',
            'carmine': '[#6]-[#6]1:[#6]2:[#6](:[#6]:[#6](:[#6]:1-[#6](=[#8])-[#8])-[#8])-[#6](=[#8])-[#6]1:[#6](-[#6]-2=[#8]):[#6](:[#6](:[#6](:[#6]:1-[#8])-[#8])-[#6]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6]-[#8])-[#8])-[#8])-[#8])-[#8]',
            'beta-carotene': '[#6]-[#6]1=[#6](-[#6](-[#6]-[#6]-[#6]-1)(-[#6])-[#6])-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6](=[#6]-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]=[#6](-[#6])-[#6]=[#6]-[#6]1=[#6](-[#6]-[#6]-[#6]-[#6]-1(-[#6])-[#6])-[#6])-[#6])-[#6]',
            'chromium hydroxide green': '[#8].[#8].[#8-2].[#8-2].[#8-2].[Cr+3].[Cr+3]',
            'chromium oxide greens': '[#8-2].[#8-2].[#8-2].[Cr+3].[Cr+3]',
            'copper, metallic powder ': '[Cu]',
            'dihydroxyacetone': '[#6](-[#6](=[#8])-[#6]-[#8])-[#8]',
            'disodium edta-copper': '[#6](-[#6]-[#7](-[#6]-[#6](=[#8])-[#8-])-[#6]-[#6](=[#8])-[#8-])-[#7](-[#6]-[#6](=[#8])-[#8-])-[#6]-[#6](=[#8])-[#8-].[Na+].[Na+].[Cu+2]',
            'ferric ammonium ferrocyanide': '[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#7H4+].[Fe+2].[Fe+3]',
            'ferric ferrocyanide': '[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[#6-]#[#7].[Fe+3].[Fe+3]',
            'guaiazulene': '[#6]-[#6]1:[#6]2:[#6]:[#6]:[#6](:[#6]-2:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#6])-[#6])-[#6]',
            'guanine': '[#6]1:[#7]:[#6]2:[#6](:[#7H]:1):[#6](=[#8]):[#7H]:[#6](:[#7]:2)-[#7]',
            'henna': '[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1)-[#6](=[#6]-[#6](=[#8])-[#6]-2=[#8])-[#8]',
            'iron oxides': '[#8]=[Fe]-[#8]-[Fe]=[#8]',
            'luminescent zinc sulfide': '[#16]=[Zn]',
            'manganese violet': '[#7H4+].[#8-]-[#15](=[#8])(-[#8-])-[#8]-[#15](=[#8])(-[#8-])-[#8-].[Mn+3]',
            'mica': '[#8-2].[#8]=[Al]-[#8]-[Al]=[#8].[#8]=[Si]=[#8].[K+].[K+]',
            'potassium sodium copper chlorophyllin': '[#6]-[#6]-[#6]1=[#6](-[#6]2:[#7]:[#6]-1:[#6]:[#6]1:[#6](:[#6](:[#6](:[#7-]:1):[#6](:[#6]1-[#6](-[#6](-[#6](:[#7]:1):[#6]:[#6]1:[#6](:[#6](:[#6](:[#6]:2):[#7-]:1)-[#6]=[#6])-[#6])-[#6])-[#6]-[#6]-[#6](=[#8])-[#8-])-[#6]-[#6](=[#8])-[#8-])-[#6](=[#8])-[#8-])-[#6])-[#6].[Na+].[Na+].[Na+].[Cu+2]',
            'pyrophyllite': '[#8]-[Si](=[#8])-[#8-].[#8]-[Si](=[#8])-[#8-].[Al+2]',
            'silver': '[Ag]',
            'silver nitrate': '[#7+](=[#8])(-[#8-])-[#8-].[Ag+]',
            'titanium dioxide': '[#8]=[Ti]=[#8]',
            'ultramarines': '[#8-]-[Si](-[#8-])(-[#8-])-[#8-].[#8-]-[Si](-[#8-])(-[#8-])-[#8-].[#8-]-[Si](-[#8-])(-[#8-])-[#8-].[#8-]-[Si](-[#8-])(-[#8-])-[#8-].[#8-]-[Si](-[#8-])(-[#8-])-[#8-].[#8-]-[Si](-[#8-])(-[#8-])-[#8-].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Al+3].[Al+3].[Al+3].[Al+3].[Al+3].[Al+3].[#16-]-[#16]-[#16-]',
            'zinc oxide': '[#8]=[Zn]',
        }

        return smarts