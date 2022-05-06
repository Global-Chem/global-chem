 #!/usr/bin/env python3
 #
 # GlobalChem - Trastuzumab Herzuma
 #
 # -----------------------------------

 class Trastuzumab_Herzuma(object):

     def __init__(self):

         self.name = 'trastuzumab_herzuma'

     @staticmethod
     def get_smiles():

         smiles = {
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'histidine hydrochloride monohydrate': 'C1=C(NC=N1)C[C@@H](C(=O)O)N.O.Cl',
             'polysorbate 20': 'CCCCCCCCCCCC(=O)OCCOCC(C1C(C(CO1)OCCO)OCCO)OCCO',
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'histidine': ''
             'histidine hydrochloride monohydrate': ''
             'polysorbate 20': '',
         }

         return smarts
