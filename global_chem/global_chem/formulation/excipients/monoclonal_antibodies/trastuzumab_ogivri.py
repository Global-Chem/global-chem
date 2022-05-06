 #!/usr/bin/env python3
 #
 # GlobalChem - Trastuzumab Ogivri
 #
 # -----------------------------------

 class Trastuzumab_Ogivri(object):

     def __init__(self):

         self.name = 'trastuzumab_ogivri'

     @staticmethod
     def get_smiles():

         smiles = {
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'histidine hydrochloride monohydrate': 'C1=C(NC=N1)C[C@@H](C(=O)O)N.O.Cl',
             'polyethylene glycol 3350/macrogol 3350': '',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'histidine': ''
             'histidine hydrochloride monohydrate': ''
             'polyethylene glycol 3350/macrogol 3350': '',
         }

         return smarts
