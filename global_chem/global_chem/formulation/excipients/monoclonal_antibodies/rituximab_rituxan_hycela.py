 #!/usr/bin/env python3
 #
 # GlobalChem - Rituximab Rituxan Hycela
 #
 # -----------------------------------

 class Rituximab_Rituxan_Hycela(object):

     def __init__(self):

         self.name = 'rituximab_rituxan_hycela'

     @staticmethod
     def get_smiles():

         smiles = {
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'histidine hydrochloride monohydrate': 'C1=C(NC=N1)C[C@@H](C(=O)O)N.O.Cl',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'trehalose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O'
             '
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'histidine': '',
             'histidine hydrochloride monohydrate': '',
             'polysorbate 80': ''
             'trehalose': ''
         }

         return smarts
