 #!/usr/bin/env python3
 #
 # GlobalChem - Basiliximab Simulect
 #
 # -----------------------------------

 class Basiliximab_Simulect(object):

     def __init__(self):

         self.name = 'basiliximab_simulect'

     @staticmethod
     def get_smiles():

         smiles = {
             'glycine': 'C(C(=O)O)N',
             'mannitol': 'C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O',
             'sodium phosphate, dibasic, anhydrous': 'OP(=O)([O-])[O-].[Na+].[Na+]',
             'potassium phosphate monobasic': 'OP(=O)(O)[O-].[K+]',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O'
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'glycine': '',
             'mannitol': '',
             'sodium phosphate, dibasic, anhydrous': '',
             'potassium phosphate monobasic': ''
             'sucrose': ''
         }

         return smarts
