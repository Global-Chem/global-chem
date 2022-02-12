#!/usr/bin/env python3
#
# GlobalChem - Interstellar Space
#
# -------------------------------

class InterstellarSpace(object):

    def __init__(self):

        pass

    @staticmethod
    def smiles():

        smiles = {
            'aluminum monochloride': '[Al]Cl',
            'aluminum monofluoride': '[Al]F',
            'aluminum isocyanide': '[Al][C-]#[NH+]',
            'methylidyne': '[CH]',
            'methyliumylidene': '[C+]',
            'hydrogen cyanide': 'CN',
            'hydrogen isocyanide': '[C-]#[NH+]',
            'isocyanic acid': 'O=C=N',
            'oxomethyl': '[C+]=O',
            'oxomethylium': '[C-]=O',
            'hydroxymethylidyne': '[C+]O',
            'hydroxyoxomethylium': 'O#[C+]O',
            'thiooxomethylium': 'S#[C+]O',
            'methylene': '[CH2]',
            'iminomethylium': 'C#[NH+]',
            'methylene amidogen': 'C=[N]',
            'cyanamide': 'C(#N)N',
            'formaldehyde': 'C=O',
            'formic acid': 'C(=O)O',
            'thioformaldehyde': 'C=S',
            'methyl': '[CH3]',
            'methanimine': 'C=N',
            'formamide': 'C(=O)N',
            'hydroxy methylium ion': 'C=[OH+]',
            'methane': 'C',
            'methanol': 'CO',
            'methanethiol': 'CS',
            'methylamine': 'CN',
            'magnesium cyanide': '[C-]#N.[C-]#N.[Mg+2]',
            'magnesium isocyanide': '[C-]#[N+].[C-]#[N+].[Mg+2]',
            'cyanide radical': '[C]#N',
            'cyanide radical ion': '[C]#[N+]',
            'sodium cyanide': '[C-]#N.[Na+]',
            'silicon cyanide': '[C-]#N.[Si+]',
            'cyanoimidogen': '[C-]#N.[N+]',
            'carbon monoxide': '[C-]#O',
            'carbon monoxide ion': '[C-]#[O+]',
            'carbon oxysulfide': '',
            'carbon dioxide': 'O=C=O',
            'carbon dioxide ion': 'O=C=[O+]',
            'carbon phosphide': '[C-]#[P+]',
            'carbon monosulfide': '[C-]#[S+]',
            'silicon carbide': '[C-]#[Si+]',
            'dicarbon': 'C#[C+]',
            'ethynyl': 'C#[C]',
            'cyanomethylene': 'CC#N',
            'acetylene': 'C#C',
            'cyanomethyl': '[CH2]C#N',
            'ketene': 'C=C=O',
            'acetonitrile': 'CC#N',
            'isocyanomethane': 'C[N+]#[C-]',
            'ethylene': 'C=C',
            'acetaldehyde': 'CC=O',
            'ethylene oxide': 'C1CO1',
            'ethenol': 'OC=C',
            'methyl formate': 'COC=O',
            'acetic acid': 'CC(=O)O',
            'glycolaldehyde': 'C(C=O)O',
            'ethane': 'CC',
            'dimethyl ether': 'COC',
            'ethylene glycol': 'C(CO)O',
            'oxoethenylidene': 'O=C=C',
            'thioxoethenylidene': 'S=C=C',
            'silicon dicarbide': '[C-]#[Si++]#[C-]',
            'tricarbon': '[C]=C=[C]',
            'cyclopropenylidyne': 'C1=CC1',
            'propenylidyne': 'C=CC',
            'cyanoacetylene': 'C#CC#N',
            'isocyanoacetylene': 'C#C[C-]#[NH+]',
            'cyclopropenylidene': 'C1=C[C]1',
            'propadienylidene': 'C=C=C',
            'protonated cyanoacetylene': 'C#CC#[NH+]',
            '2-propynal': 'C#CC=O',
            'acrylonitrile': 'C=CC#N',
            'propyne': 'CC#C',
            'propanenitrile': 'N#CCC',
            'acetone': 'CC(=O)C',
            'cyanoethynl': 'C[CH]C#N',
            '1,2-propadienylidene, 3-oxo': 'C=CCO',
            '1,2-propadienylidene, thioxo': 'C=CCS',
            'silicon tricarbon': '[C]=C=C.[Si+]',
            '1,3-butadiynyl': 'C#CC#C',
            'butatrienylidene': 'C=C=C=[C]',
            '2-butynenitrile': 'CC#CC#N',
            'silicon tetracarbide': '[C]=C=C=C.[Si+]',
            'pentacarbon': 'C=C=C=C=C',
            '2,4-pentadiynylidyne': '',
            '2,4-pentadiynenitrile': '',
            '1,3-pentadiyne': '',
            '1,3-butadiynylium, 4-cyano': '',
            '1,3,5-hexatriynyl': '',
            '1,3,5-hexatriyne': '',
            '1,2,3,4,5-hexapentaenylidene': '',
            'benzene': 'C1=CC=CC=C1',
            '2,4,6-heptatriynylidyne': '',
            '2,4,6-heptatriynenitrile': '',
            '1,3,5,7-octatetraynyl': '',
            '2,4,6,8-nonatetraynenitrile': '',
            '2,4,6,8,10-undecapentaynenitrile': '',
            'hydrogen chloride': 'HCl',
            'potassium chloride': '[Cl-].[K+]',
            'sodium chloride': '[Cl-].[Na+]',
            'hydrogen fluoride': 'HF',
            'iron monoxide': 'O.[Fe]',
            'lithium hydride': '[H-].[Li+]',
            'imidogen': '[NH]',
            'nitrosyl hydride': 'N=O',
            'hydrodinitrogen': 'N#[NH+]',
            'hydroxyl': 'OH',
            'oxoniumylidene': '',
            'mercapto': '[SH]',
            'hydrogen': '[H]',
            'amidogen': '[NH2]',
            'water': 'O',
            'oxoniumyl': '[OH2+]',
            'hydrogen sulfide': 'S',
            'trihydrogen ion': '[H+]1[H][H]1',
            'ammonia': 'N',
            'oxonium hydride': '[OH3+]',
            'silane': '[SiH4]',
            'nitric oxide': '[N]=O',
            'phosphorous nitride': 'N#P',
            'nitrogen sulfide': '[N]=S',
            'silicon nitride': 'N12[Si]34N5[Si]16N3[Si]25N46',
            'nitrogen ion': 'N#N',
            'nitrous oxide': '[N-]=[N+]=O',
            'sulfur monoxide': 'O=S',
            'silicon monoxide': 'O.[Si]',
            'sulfur dioxide': 'O=S=O',
            'silicon monosulfide': '[Si]=S',
            'disulfur': 'S=S'
        }

        return smiles

    @staticmethod
    def get_smarts():


        smarts = {
            'aluminum monochloride': '[Al]-[#17]',
            'aluminum monofluoride': '[Al]-[#9]',
            'methylidyne': '[#6H]',
            'methyliumylidene': '[#6+]',
            'hydrogen cyanide': '[#6]-[#7]',
            'hydrogen isocyanide': '[#6-]#[#7H+]',
            'isocyanic acid': '[#8]=[#6]=[#7]',
            'oxomethyl': '[#6+]=[#8]',
            'oxomethylium': '[#6-]=[#8]',
            'hydroxymethylidyne': '[#6+]-[#8]',
            'thiooxomethylium': '[#16]#[#6+]-[#8]',
            'methylene': '[#6H2]',
            'iminomethylium': '[#6]#[#7H+]',
            'methylene amidogen': '[#6]=[#7]',
            'cyanamide': '[#6](#[#7])-[#7]',
            'formaldehyde': '[#6]=[#8]',
            'formic acid': '[#6](=[#8])-[#8]',
            'thioformaldehyde': '[#6]=[#16]',
            'methyl': '[#6H3]',
            'methanimine': '[#6]=[#7]',
            'formamide': '[#6](=[#8])-[#7]',
            'hydroxy methylium ion': '[#6]=[#8H+]',
            'methane': '[#6]',
            'methanol': '[#6]-[#8]',
            'methanethiol': '[#6]-[#16]',
            'methylamine': '[#6]-[#7]',
            'magnesium cyanide': '[#6-]#[#7].[#6-]#[#7].[Mg+2]',
            'magnesium isocyanide': '[#6-]#[#7+].[#6-]#[#7+].[Mg+2]',
            'cyanide radical': '[#6]#[#7]',
            'cyanide radical ion': '[#6]#[#7+]',
            'sodium cyanide': '[#6-]#[#7].[Na+]',
            'silicon cyanide': '[#6-]#[#7].[Si+]',
            'cyanoimidogen': '[#6-]#[#7].[#7+]',
            'carbon monoxide ion': '[#6-]#[#8+]',
            'carbon oxysulfide': '',
            'carbon dioxide': '[#8]=[#6]=[#8]',
            'carbon dioxide ion': '[#8]=[#6]=[#8+]',
            'carbon phosphide': '[#6-]#[#15+]',
            'carbon monosulfide': '[#6-]#[#16+]',
            'silicon carbide': '[#6-]#[Si+]',
            'dicarbon': '[#6]#[#6+]',
            'ethynyl': '[#6]#[#6]',
            'cyanomethylene': '[#6]-[#6]#[#7]',
            'acetylene': '[#6]#[#6]',
            'cyanomethyl': '[#6H2]-[#6]#[#7]',
            'ketene': '[#6]=[#6]=[#8]',
            'acetonitrile': '[#6]-[#6]#[#7]',
            'isocyanomethane': '[#6]-[#7+]#[#6-]',
            'ethylene': '[#6]=[#6]',
            'acetaldehyde': '[#6]-[#6]=[#8]',
            'ethylene oxide': '[#6]1-[#6]-[#8]-1',
            'ethenol': '[#8]-[#6]=[#6]',
            'methyl formate': '[#6]-[#8]-[#6]=[#8]',
            'acetic acid': '[#6]-[#6](=[#8])-[#8]',
            'glycolaldehyde': '[#6](-[#6]=[#8])-[#8]',
            'ethane': '[#6]-[#6]',
            'dimethyl ether': '[#6]-[#8]-[#6]',
            'ethylene glycol': '[#6](-[#6]-[#8])-[#8]',
            'oxoethenylidene': '[#8]=[#6]=[#6]',
            'thioxoethenylidene': '[#16]=[#6]=[#6]',
            'silicon dicarbide': '[#6-]#[Si+2]#[#6-]',
            'tricarbon': '[#6]=[#6]=[#6]',
            'cyclopropenylidyne': '[#6]1=[#6]-[#6]-1',
            'propenylidyne': '[#6]=[#6]-[#6]',
            'cyanoacetylene': '[#6]#[#6]-[#6]#[#7]',
            'cyclopropenylidene': '[#6]1=[#6]-[#6]-1',
            'propadienylidene': '[#6]=[#6]=[#6]',
            'protonated cyanoacetylene': '[#6]#[#6]-[#6]#[#7H+]',
            '2-propynal': '[#6]#[#6]-[#6]=[#8]',
            'acrylonitrile': '[#6]=[#6]-[#6]#[#7]',
            'propyne': '[#6]-[#6]#[#6]',
            'propanenitrile': '[#7]#[#6]-[#6]-[#6]',
            'acetone': '[#6]-[#6](=[#8])-[#6]',
            'cyanoethynl': '[#6]-[#6H]-[#6]#[#7]',
            '1,2-propadienylidene, 3-oxo': '[#6]=[#6]-[#6]-[#8]',
            '1,2-propadienylidene, thioxo': '[#6]=[#6]-[#6]-[#16]',
            'silicon tricarbon': '[#6]=[#6]=[#6].[Si+]',
            '1,3-butadiynyl': '[#6]#[#6]-[#6]#[#6]',
            'butatrienylidene': '[#6]=[#6]=[#6]=[#6]',
            '2-butynenitrile': '[#6]-[#6]#[#6]-[#6]#[#7]',
            'silicon tetracarbide': '[#6]=[#6]=[#6]=[#6].[Si+]',
            'pentacarbon': '[#6]=[#6]=[#6]=[#6]=[#6]',
            '2,4-pentadiynylidyne': '',
            '2,4-pentadiynenitrile': '',
            '1,3-pentadiyne': '',
            '1,3-butadiynylium, 4-cyano': '',
            '1,3,5-hexatriynyl': '',
            '1,3,5-hexatriyne': '',
            '1,2,3,4,5-hexapentaenylidene': '',
            'benzene': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
            '2,4,6-heptatriynylidyne': '',
            '2,4,6-heptatriynenitrile': '',
            '1,3,5,7-octatetraynyl': '',
            '2,4,6,8-nonatetraynenitrile': '',
            '2,4,6,8,10-undecapentaynenitrile': '',
            'potassium chloride': '[#17-].[K+]',
            'sodium chloride': '[#17-].[Na+]',
            'iron monoxide': '[#8].[Fe]',
            'lithium hydride': '[H-].[Li+]',
            'imidogen': '[#7H]',
            'nitrosyl hydride': '[#7]=[#8]',
            'hydrodinitrogen': '[#7]#[#7H+]',
            'oxoniumylidene': '',
            'mercapto': '[#16H]',
            'hydrogen': '[H]',
            'amidogen': '[#7H2]',
            'water': '[#8]',
            'oxoniumyl': '[#8H2+]',
            'hydrogen sulfide': '[#16]',
            'ammonia': '[#7]',
            'oxonium hydride': '[#8H3+]',
            'silane': '[SiH4]',
            'nitric oxide': '[#7]=[#8]',
            'phosphorous nitride': '[#7]#[#15]',
            'nitrogen sulfide': '[#7]=[#16]',
            'silicon nitride': '[#7]12-[Si]34-[#7]5-[Si]-16-[#7]-3-[Si]-2-5-[#7]-4-6',
            'nitrogen ion': '[#7]#[#7]',
            'nitrous oxide': '[#7-]=[#7+]=[#8]',
            'sulfur monoxide': '[#8]=[#16]',
            'silicon monoxide': '[#8].[Si]',
            'sulfur dioxide': '[#8]=[#16]=[#8]',
            'silicon monosulfide': '[Si]=[#16]',
            'disulfur': '[#16]=[#16]',
        }

        return smarts