 #!/usr/bin/env python3
 #
 # GlobalChem - Palivizumab Synagis
 #
 # -----------------------------------

 class Palivizumab_Synagis(object):

     def __init__(self):

         self.name = 'palivizumab_synagis'

     @staticmethod
     def get_smiles():

         smiles = {
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'glycine': 'C(C(=O)O)N',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'histidine': '',
             'glycine': '',
         }

         return smarts
