 #!/usr/bin/env python3
 #
 # GlobalChem - Trastuzumab Enhertu
 #
 # -----------------------------------

 class Trastuzumab_Enhertu(object):

     def __init__(self):

         self.name = 'trastuzumab_enhertu'

     @staticmethod
     def get_smiles():

         smiles = {
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'histidine hydrochloride monohydrate': 'C1=C(NC=N1)C[C@@H](C(=O)O)N.O.Cl',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
             'd-sorbitol': 'C([C@H]([C@H]([C@@H]([C@H](CO)O)O)O)O)O'
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'histidine': ''
             'histidine hydrochloride monohydrate': ''
             'polysorbate 80': '',
             'sucrose': '',
             'd-sorbitol': ''
         }

         return smarts
