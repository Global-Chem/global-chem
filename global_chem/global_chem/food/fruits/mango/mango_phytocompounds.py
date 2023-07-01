#!/usr/bin/env python3
#
# GlobalChem - Common Salts
#
# -----------------------------------

class Mango(object):


    '''

    Missing Entries:

          'alpha phellandrene': '',
          'alpha humulene': '',
          'gamma terpinene': '',
          'alpha pinene': '',
          'trans caryophyllene': '',
          'sabinene': '',
          'carene': '',
          'cis-caryophyllene': '',
          'αlpha humulene': '',
          'germacrene d': '',
          'aromadendrene': '',
          'beta ubebene': '',
          'alpha cubebene': '',
          'alpha bourbonene': '',
          'beta elemene': '',
          'hepta-2,4(E,E)-dienoic acid': '',
          'alpha-Linoleic acid': '',
    '''

    def __init__(self):

      self.name = 'mango'

    @staticmethod
    def get_smiles():

        smiles = {
          'arginine': 'C(CC(C(=O)O)N)CN=C(N)N',
          'cysteine': 'C([C@@H](C(=O)O)N)S',
          'glycine': 'C(C(=O)O)N',
          'isoleucine': 'CC[C@H](C)[C@@H](C(=O)O)N ',
          'leucine': 'CC(C)C[C@@H](C(=O)O)N',
          'lysine': 'C(CCN)C[C@@H](C(=O)O)N',
          'methionine': 'CSCC[C@@H](C(=O)O)N',
          'phenylalanine': 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N',
          'proline': 'C1C[C@H](NC1)C(=O)O',
          'serine': 'C([C@@H](C(=O)O)N)O',
          'threonine': 'C[C@H]([C@@H](C(=O)O)N)O',
          'tryptophan': 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N ',
          'tyrosine': 'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O  ',
          'valine': 'CC(C)[C@@H](C(=O)O)N',
          'glutamic acid': 'C(CC(=O)O)[C@@H](C(=O)O)N',
          'aspartic acid': 'C([C@@H](C(=O)O)N)C(=O)O',
          'alanine': 'C[C@@H](C(=O)O)N',
          'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
          'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
          'stearic acid ': 'CCCCCCCCCCCCCCCCCC(=O)O',
          'lignoceric acid': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)O',
          'oleic acid': r'CCCCCCCC/C=C\CCCCCCCC(=O)O',
          'linoleic acid': r'CCCCC/C=C\C/C=C\CCCCCCCC(=O)O',
          'myristic acid': 'CCCCCCCCCCCCCC(=O)O',
          'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O ',
          'arachidic acid': 'CCCCCCCCCCCCCCCCCCCC(=O)O',
          'behenic acid': 'CCCCCCCCCCCCCCCCCCCCCC(=O)O',
          'palmitoleic acid': r'CCCCCC/C=C\CCCCCCCC(=O)O',
          'hexadecenoic acid': r'CCCC/C=C/CCCCCCCCCC(=O)O',
          'heptadecenoic acid': r'CCCCCC/C=C/CCCCCCCCC(=O)O',
          'octadecenoic acid': r'CCCCCC/C=C/CCCCCCCCCC(=O)O',
          'eicosenoic acid': r'CCCCCCCC/C=C\CCCCCCCCCC(=O)O',
          '9,12-Hexadecadienoic acid': r'CCC/C=C/C/C=C/CCCCCCCC(=O)O',
          '9,15-Octadecadienoic acid': r'CC/C=C/CCCC/C=C/CCCCCCCC(=O)O',
          'linolenic acid': r'CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O',
          'ascorbic acid ': 'C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O',
          'thiamine': 'CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO',
          'riboflavin': 'CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@@H]([C@@H]([C@@H](CO)O)O)O',
          'niacin': 'C1=CC(=CN=C1)C(=O)O',
          'pantothenic acid': 'CC(C)(CO)[C@H](C(=O)NCCC(=O)O)O',
          'pyridoxine ': 'CC1=NC=C(C(=C1O)CO)CO',
          'folic acid ': 'C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N',
          'vitamin A': r'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C',
          'vitamin E': 'CC1=C(C2=C(CC[C@@](O2)(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=C1O)C)C',
          'vitamin K': 'CC1=C(C(=O)C2=CC=CC=C2C1=O)C/C=C(\C)/CCCC(C)CCCC(C)CCCC(C)C',
          'gallic Acid': 'C1=C(C=C(C(=C1O)O)O)C(=O)O',
          'vanillic acid ': 'COC1=C(C=CC(=C1)C(=O)O)O',
          'syringic acid ': 'COC1=CC(=CC(=C1O)OC)C(=O)O',
          'protocatechuic acid': 'C1=CC(=C(C=C1C(=O)O)O)O',
          'para hydroxybenzoic acid ': 'C1=CC(=CC=C1C(=O)O)O',
          'paracoumaric acid': r'C1=CC(=CC=C1/C=C/C(=O)O)O',
          'chlorogenic acid ': r'C1[C@H]([C@H]([C@@H](C[C@@]1(C(=O)O)O)OC(=O)/C=C/C2=CC(=C(C=C2)O)O)O)O',
          'ferulic acid': r'COC1=C(C=CC(=C1)/C=C/C(=O)O)O',
          'caffeic acid': r'C1=CC(=C(C=C1/C=C/C(=O)O)O)O',
          'theogallin': 'C1[C@H]([C@H]([C@@H](C[C@@]1(C(=O)O)O)OC(=O)C2=CC(=C(C(=C2)O)O)O)O)O',
          'quercetin-3-O-galactoside': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O)O)O',
          'quercetin-3-O-glucoside': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O',
          'quercetin-3-O-xyloside': 'C1C(C(C(C(O1)OC2=C(OC3=CC(=CC(=C3C2=O)O)O)C4=CC(=C(C=C4)O)O)O)O)O',
          'magniferin': 'C1=C2C(=CC(=C1O)O)OC3=C(C2=O)C(=C(C(=C3)O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O',
          'cyanidin': 'C1=CC(=C(C=C1C2=[O+]C3=CC(=CC(=C3C=C2O)O)O)O)O',
          'delphinidin': 'C1=C(C=C(C(=C1O)O)O)C2=[O+]C3=CC(=CC(=C3C=C2O)O)O.[Cl-]',
          'pelargonidin': 'C1=CC(=CC=C1C2=[O+]C3=CC(=CC(=C3C=C2O)O)O)O',
          'catechin': 'C1[C@@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O',
          'apigenin': 'C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O ',
          'luteolin': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O',
          'kaempferol': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O',
          'myricetin': 'C1=C(C=C(C(=C1O)O)O)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O',
          '9-cis-violaxanthin': r'C/C(=C\C=C\C=C(/C)\C=C\C=C(\C)/C=C/[C@]12[C@](O1)(C[C@H](CC2(C)C)O)C)/C=C/C=C(\C)/C=C/[C@]34[C@](O3)(C[C@H](CC4(C)C)O)C',
          'limonene': 'CC1=CCC(CC1)C(=C)C',
          'alpha-terpinolene': r'C/C(C)=C1CCC(C)C=C\1',
          'd-carvone': 'CC1=CC[C@@H](CC1=O)C(=C)C',
        }

        return smiles


