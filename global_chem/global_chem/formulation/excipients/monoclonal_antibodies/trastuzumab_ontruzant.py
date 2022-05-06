 #!/usr/bin/env python3
 #
 # GlobalChem - Trastuzumab Ontruzant
 #
 # -----------------------------------

 class Trastuzumab_Ontruzant(object):

     def __init__(self):

         self.name = 'trastuzumab_ontruzant'

     @staticmethod
     def get_smiles():

         smiles = {
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'histidine hydrochloride monohydrate': 'C1=C(NC=N1)C[C@@H](C(=O)O)N.O.Cl',
             'polysorbate 20': 'CCCCCCCCCCCC(=O)OCCOCC(C1C(C(CO1)OCCO)OCCO)OCCO',
             'trehalose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'histidine': ''
             'histidine hydrochloride monohydrate': ''
             'polysorbate 20': '',
             'trehalose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O',
         }

         return smarts
