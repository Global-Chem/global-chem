 #!/usr/bin/env python3
 #
 # GlobalChem - Basiliximab Simulect
 #
 # -----------------------------------

 class BasiliximabSimulect(object):

     def __init__(self):

         self.name = 'basiliximab_simulect'

     @staticmethod
     def get_smiles():

         smiles = {
             'glycine': 'C(C(=O)O)N',
             'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
             'sodium phosphate, dibasic, anhydrous': 'OP(=O)([O-])[O-].[Na+].[Na+]',
             'potassium phosphate monobasic': 'OP(=O)(O)[O-].[K+]',
             'sucrose': 'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'glycine': '',
             'mannitol': '',
             'sodium phosphate, dibasic, anhydrous': '',
             'potassium phosphate monobasic': '',
             'sucrose': '',
         }

         return smarts
