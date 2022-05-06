#!/usr/bin/env python3
 #
 # GlobalChem - Rituximab Ruxience
 #
 # -----------------------------------

 class Rituximab_Ruxience(object):

     def __init__(self):

         self.name = 'rituximab_ruxience'

     @staticmethod
     def get_smiles():

         smiles = {
             'edetate disodium dihydrate': 'C(CN(CC(=O)O)CC(=O)[O-])N(CC(=O)O)CC(=O)[O-].O.O.[Na+].[Na+]',
             'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
             'histidine hydrochloride monohydrate': 'C1=C(NC=N1)C[C@@H](C(=O)O)N.O.Cl',
             'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
             'sucrose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O'
             '
         }

         return smiles

     @staticmethod
     def get_smarts():

         smarts = {
             'edetate disodium dihydrate': '',
             'histidine': '',
             'histidine hydrochloride monohydrate': '',
             'polysorbate 80': ''
             'sucrose': ''
         }

         return smarts
