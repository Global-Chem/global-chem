#!/usr/bin/env python3
#
# GlobalChem - Ginger Terpenoids
# ------------------------------

class ThaiGingerTerpenoids(object):

  def __init__(self):

      self.name = 'thai_ginger_terpenoids'

  @staticmethod
  def get_smiles():

      '''

      Missing:

          '6β, 14β-dihydroxyisopimara-8(9),15-diene' : '',

      '''

      smiles = {
        '3-caren-5-one'                                      : 'CC1=CC(=O)C2C(C1)C2(C)C',
        '(3r,4r,6s)-3,6-dihydroxy-1-menthene'                : 'CC1=C[C@@H]([C@H](C[C@@H]1O)C(C)C)O',
        '(1R,2S,4R)-p-menth-5-ene-1,2,8-triol'               : 'O[C@H]1[C@@H]([C@H](C=C)[C@H](O)C1)CO',
        'oxyphyllenodiol B'                                  : 'CC(C)[C@H]1CCC(=O)C2=C1[C@@H]([C@](CC2)(C)O)O',
        'hedytriol'                                          : 'CC(C)(C)C1=CC(=C(C(=C1)O)O)C(C)(C)C',
        'kaemgalangol A'                                     : 'CC(C)(C)CCCC(C)(C)CCCC(C)(C)CCCC(C)(C)CCCC1=CC(=C(C(=C1)O)O)C(C)(C)C',
        '6β-hydroxypimara-8(14), 15-diene-1-one'             : 'O[C@@H]1[C@@H](C)C[C@@H]2[C@]1(CC[C@@H]3[C@H]2CC=C4[C@@]3(C)CCC(=C4)C(=O)C)C',
        'sandaracopimaradien-6β, 9α-diol-l-one'              : 'O=C[C@@]1(C)[C@H]2[C@@H](O)C[C@H]3[C@](O)(C)[C@]4(C)CC[C@@H]5[C@@]4(C)CC[C@@H](O)C(C)(C)[C@]35C2=C1',
        '(-)-sandaracopimaradiene'                           : 'C[C@@]1(CC[C@H]2C(=C1)CC[C@@H]3[C@@]2(CCCC3(C)C)C)C=C',
        'sandaracopimaradiene-9α-ol'                         : 'CC(C)[C@]1(C)CC[C@@H](CC1)C(C)C2=CC(=C[C@@H](O)C2)C',
        'kaempulchraol l'                                    : 'C[C@@]1(CC[C@]2(C(=C1)C[C@H]([C@@H]3[C@@]2(CCCC3(C)C)C)O)OC)C=C',
        'kaempulchraol E'                                    : 'C[C@]12[C@H]3CC[C@@H](C=C3C[C@H]([C@@H]1C(CC[C@@H]2O)(C)C)O)C=C',
        '8(14),15-sandaracopimaradiene- 1a,9α-diol'          : 'CC1=C2[C@]([C@H]3[C@]4([C@@H]([C@H](C)C(C)(C)[C@H]4CC[C@@]3(C)C2)CC[C@H]1O)C)(C)C',
        'kaempulchraol L'                                    : 'C[C@@]1(CC[C@]2(C(=C1)C[C@H]([C@@H]3[C@@]2(CCCC3(C)C)C)O)OC)C=', # Mistake
        '2α-acetoxy sandaracopimaradien-1α-ol'               : 'CC(=O)OC1(C2[C@@]([C@@H]3[C@]4([C@H](C(C)=C)[C@@H]([C@H]4CC[C@]3(C)C2)CC[C@@H]1O)(C)C)C',
        '1,11-dihydroxypimara-8(14), 15-diene'               : 'CC1=CCCC2C(C1)CCC(C2)(O)C(C=C)C(CO)C',
        '6β, 14α-dihydroxyisopimara-8(9),15-diene'           : 'CC(C)[C@]1(C)[C@@H]2CC=C3[C@@H](O)C[C@@H](O)C4=C(C)[C@H](O)C[C@@]4(C)[C@]3(C)CC[C@@H]2[C@@]1(C)C',
        '6α, 14β-dihydroxyisopimara-8(9), 15-diene'          : 'CC(C)[C@]1(C)[C@@H]2CC=C3[C@H](O)C[C@H](O)[C@@]4(C)[C@H](O)C[C@@]4(C)[C@]3(C)CC[C@@H]2[C@@]1(C)C',
        '1α-hydroxy-14α-methoxyisopimara-8(9),15-diene'      : 'CC(C)[C@]1(C)[C@@H]2CC=C3[C@H](O)[C@H](OC)[C@@]4(C)[C@H](C)C[C@@]4(C)[C@]3(C)CC[C@@H]2[C@@]1(C)C',
        '1α,14α-dihydroxyisopimara-8(9),15-diene'            : 'CC(C)[C@]1(C)[C@@H]2CC=C3[C@H](O)[C@@H](O)[C@@]4(C)[C@H](C)C[C@@]4(C)[C@]3(C)CC[C@@H]2[C@@]1(C)C',
        'boesenberol l'                                      : 'CC(C)[C@]1(C)[C@@H]2CC=C3[C@H](C)[C@@H](O)[C@@]4(C)[C@H](C)C[C@@]4(C)[C@]3(C)CC[C@@H]2[C@@]1(C)C',
        'boesenberol J'                                      : 'CC(C)[C@]1(C)[C@@H]2CC=C3[C@H](C)[C@@H](O)[C@@]4(C)[C@H](C)C[C@@]4(C)[C@]3(C)CC[C@@H]2[C@@]1(C)C',
        '6β-acetoxysandaracopimaradiene-9α-ol'               : 'CC(=O)O[C@@]1(C)[C@H]2CC[C@](C)(O)C(C)(C)[C@@]2(C)[C@@H]3CC[C@@]4(C)[C@H](O)CC[C@]4(C)[C@@]3(C)CC1',
        '6β-acetoxysandaracopimaradiene-9α-ol-1-one'         : 'CC(=O)[C@@]1(C)[C@H]2CC[C@](C)(O)C(C)(C)[C@@]2(C)[C@@H]3CC[C@@]4(C)[C@H](O)CC[C@]4(C)[C@@]3(C)CC1',
        '6β-acetoxysandaracopimaradiene-1α,9α-diol'          : 'CC(=O)[C@]1(C)[C@@H]2CC[C@](C)(O)C(C)(C)[C@@]2(C)[C@@H]3CC[C@@]4(C)[C@@H](O)CC[C@]4(C)[C@@]3(C)CC1',
        '6β-acetoxy-1α-14α-dihydroxyisopimara-8(9),15-diene' : 'CC(=O)[C@]1(C)[C@@H]2CC[C@](C)(O)[C@@]3(C)[C@H](O)C[C@@]4(C)[C@H](O)C[C@@]4(C)[C@@]3(C)CC[C@@H]2[C@@]1(C)C',
      }

      return smiles

  @staticmethod
  def get_smarts():

      smarts = {

      }

      return smarts
