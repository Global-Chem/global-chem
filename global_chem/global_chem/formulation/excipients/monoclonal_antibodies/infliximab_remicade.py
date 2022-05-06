 #!/usr/bin/env python3
 #
 # GlobalChem - Infliximab Remicade
 #
 # -----------------------------------

 class Infliximab_Remicade(object):

     def __init__(self):

         self.name = 'infliximab_remicade'

     @staticmethod
     def get_smiles():

         smiles = {
             'sodium phosphate, monobasic monohydrate': '[O.OP(=O)(O)[O-].[Na+]',
             'sodium phosphate, dibasic, dihydrate': 'O.O.OP(=O)([O-])[O-].[Na+].[Na+]',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'sodium phosphate, monobasic monohydrate': '',
             'sodium phosphate, dibasic, dihydrate': '',
             'polysorbate 80': '',
             'sucrose': '',
         }

         return smarts
