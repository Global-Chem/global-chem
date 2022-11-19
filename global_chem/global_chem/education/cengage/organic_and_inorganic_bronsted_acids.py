#!/usr/bin/env python3
#
# GlobalChem - Organic and Inorganic Bronsted Acids
#
# -------------------------------------------------

class OrganicAndInorganicBronstedAcids(object):

    def __init__(self):

        self.name = 'organic_and_inorganic_bronsted_acids'

    @staticmethod
    def get_smiles():

        smiles = {
            'sulfuric acid': 'OS(O)(=O)=O',
            'hydrogen iodide': '[H]I',
            'hydrogen bromide': '[H]Br',
            'hydrogen chloride': '[H]Cl',
            'benzenesulfonic acid': 'OS(C1=CC=CC=C1)(=O)=O',
            'methanesulfonic acid': 'CS(O)(=O)=O',
            'nitric acid': '[O-][N+](O)=O',
            'trifluoroacetic acid': 'OC(C(F)(F)F)=O',
            'benzoic acid': 'OC(C1=CC=CC=C1)=O',
            'acetic acid': 'CC(O)=O',
            'carbonic acid': 'OC(O)=O',
            'thiophenol': 'SC1=CC=CC=C1',
            'hydrogen sulfide': '[H]S[H]',
            'peracetic acid': 'CC(OO)=O',
            'phthalimide': 'O=C1NC(C2=CC=CC=C21)=O',
            'nitroethane': 'CC[N+]([O-])=O',
            'pentane-2,4-dione': 'CC(CC(C)=O)=O',
            'hydrogen cyanide': 'N#C[H]',
            'hexafluoroisopropanol': 'OC(C(F)(F)F)C(F)(F)F',
            'phenol': 'OC1=CC=CC=C1',
            'methanethiol': 'CS',
            'diethyl malonate': 'O=C(CC(OC)=O)OC',
            'cyclopentadiene': 'C1=CC=CC1',
            'water': 'O',
            'ethanol': 'OCC',
            'cyclohexanone': 'O=C1CCCCC1',
            'acetamide': 'CC(N)=O',
            'isopropanol': 'CC(C)O',
            't-butanol': 'CC(C)(C)O',
            'acetone': 'CC(C)=O',
            'ethyl acetate': 'CC(OCC)=O',
            'ethyne': 'C#C',
            'acetonitrile': 'N#CC',
            'dimethylsulfone': 'CS(C)(=O)=O',
            'dimethylsulfoxide': 'CS(C)=O',
            'ammonia': 'N',
            'diisopropylamine': 'CC(C)NC(C)C',
            'toluene': 'CC1=CC=CC=C1',
            'benzene': 'C1=CC=CC=C1',
            'propene': 'C=CC',
            'ethene': 'C=C',
            'methane': 'C',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'sulfuric acid': '[#8]-[#16](-[#8])(=[#8])=[#8]',
            'hydrogen iodide': '[#53H]',
            'hydrogen bromide': '[#35H]',
            'hydrogen chloride': '[#17H]',
            'benzenesulfonic acid': '[#8]-[#16](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(=[#8])=[#8]',
            'methanesulfonic acid': '[#6]-[#16](-[#8])(=[#8])=[#8]',
            'nitric acid': '[#8-]-[#7+](-[#8])=[#8]',
            'trifluoroacetic acid': '[#8]-[#6](-[#6](-[#9])(-[#9])-[#9])=[#8]',
            'benzoic acid': '[#8]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8]',
            'acetic acid': '[#6]-[#6](-[#8])=[#8]',
            'carbonic acid': '[#8]-[#6](-[#8])=[#8]',
            'thiophenol': '[#16]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'hydrogen sulfide': '[#16H2]',
            'peracetic acid': '[#6]-[#6](-[#8]-[#8])=[#8]',
            'phthalimide': '[#8]=[#6]1-[#7]-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1)=[#8]',
            'nitroethane': '[#6]-[#6]-[#7+](-[#8-])=[#8]',
            'pentane-2,4-dione': '[#6]-[#6](-[#6]-[#6](-[#6])=[#8])=[#8]',
            'hydrogen cyanide': '[#7]#[#6H]',
            'hexafluoroisopropanol': '[#8]-[#6](-[#6](-[#9])(-[#9])-[#9])-[#6](-[#9])(-[#9])-[#9]',
            'phenol': '[#8]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'methanethiol': '[#6]-[#16]',
            'diethyl malonate': '[#8]=[#6](-[#6]-[#6](-[#8]-[#6])=[#8])-[#8]-[#6]',
            'cyclopentadiene': '[#6]1=[#6]-[#6]=[#6]-[#6]-1',
            'water': '[#8]',
            'ethanol': '[#8]-[#6]-[#6]',
            'cyclohexanone': '[#8]=[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
            'acetamide': '[#6]-[#6](-[#7])=[#8]',
            'isopropanol': '[#6]-[#6](-[#6])-[#8]',
            't-butanol': '[#6]-[#6](-[#6])(-[#6])-[#8]',
            'acetone': '[#6]-[#6](-[#6])=[#8]',
            'ethyl acetate': '[#6]-[#6](-[#8]-[#6]-[#6])=[#8]',
            'ethyne': '[#6]#[#6]',
            'acetonitrile': '[#7]#[#6]-[#6]',
            'dimethylsulfone': '[#6]-[#16](-[#6])(=[#8])=[#8]',
            'dimethylsulfoxide': '[#6]-[#16](-[#6])=[#8]',
            'ammonia': '[#7]',
            'diisopropylamine': '[#6]-[#6](-[#6])-[#7]-[#6](-[#6])-[#6]',
            'toluene': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'benzene': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            'propene': '[#6]=[#6]-[#6]',
            'ethene': '[#6]=[#6]',
            'methane': '[#6]',
        }

        return smarts
