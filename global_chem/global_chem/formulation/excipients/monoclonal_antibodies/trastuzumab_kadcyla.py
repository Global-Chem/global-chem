 #!/usr/bin/env python3
 #
 # GlobalChem - Trastuzumab Kadcyla
 #
 # -----------------------------------

 class Trastuzumab_Kadcyla(object):

     def __init__(self):

         self.name = 'trastuzumab_kadcyla'

     @staticmethod
     def get_smiles():

         smiles = {
             'polysorbate 20': 'CCCCCCCCCCCC(=O)OCCOCC(C1C(C(CO1)OCCO)OCCO)OCCO',
             'succinic acid': 'C(CC(=O)O)C(=O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'polysorbate 20': '',
             'succinic acid': '',
             'sucrose': '',
         }

         return smarts
