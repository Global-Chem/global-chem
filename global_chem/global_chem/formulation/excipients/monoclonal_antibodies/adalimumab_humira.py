 #!/usr/bin/env python3
 #
 # GlobalChem - Adalimumab Humira
 #
 # -----------------------------------

 class Adalimumab_Humira(object):

     def __init__(self):

         self.name = 'adalimumab_humira'

     @staticmethod
     def get_smiles():

         smiles = {
             'mannitol': 'C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O'
             '
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'mannitol': '',
             'polysorbate 80': ''
             'sucrose': ''
         }

         return smarts
