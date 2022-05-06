 #!/usr/bin/env python3
 #
 # GlobalChem - Infliximab Avsola
 #
 # -----------------------------------

 class Infliximab_Avsola(object):

     def __init__(self):

         self.name = 'infliximab_avsola'

     @staticmethod
     def get_smiles():

         smiles = {
             'sodium phosphate, dibasic, anhydrous': 'OP(=O)([O-])[O-].[Na+].[Na+]',
             'sodium phosphate, monobasic, monohydrate': 'OP(=O)(O)[O-].[Na+]',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'sodium phosphate, dibasic, anhydrous': '',
             'sodium phosphate, monobasic, monohydrate': '',
             'polysorbate 80': '',
             'sucrose': '',
         }

         return smarts
