 #!/usr/bin/env python3
 #
 # GlobalChem - Infliximab Infliximab
 #
 # -----------------------------------

 class Infliximab_Infliximab(object):

     def __init__(self):

         self.name = 'infliximab_infliximab'

     @staticmethod
     def get_smiles():

         smiles = {
             'sodium phosphate, dibasic, dihydrate': 'O.O.OP(=O)([O-])[O-].[Na+].[Na+]',
             'sodium phosphate, monobasic, monohydrate': '[O.OP(=O)(O)[O-].[Na+]',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'sodium phosphate, dibasic, dihydrate': '',
             'sodium phosphate, monobasic, monohydrate': '',
             'polysorbate 80': '',
             'sucrose': '',
         }

         return smarts
