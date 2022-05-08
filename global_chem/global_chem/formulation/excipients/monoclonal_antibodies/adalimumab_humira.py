 #!/usr/bin/env python3
 #
 # GlobalChem - Adalimumab Humira
 #
 # -----------------------------------

 class AdalimumabHumira(object):

     def __init__(self):

         self.name = 'adalimumab_humira'

     @staticmethod
     def get_smiles():

         smiles = {
             'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
             'polysorbate 80': r'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O'
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
