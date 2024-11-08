#!/usr/bin/env python3
#
# GlobalChem - Benzimidazole Opioids
# Reference: https://en.wikipedia.org/wiki/List_of_benzimidazole_opioids
# ------------------------------------------

class BenzimidazoleOpioids(object):
     
    def __init__(self):

        self.name = 'benzimidazole_opioids'

    @staticmethod
    def get_smiles():

         smiles = {
         'Desnitazene (1-diethylaminoethyl-2-benzyl-benzimidazole)' : 'CCN(CC)CCN1C(CC2=CC=CC=C2)=NC3=CC=CC=C31',
         'Metodesnitazene(Metazene)' : 'CCN(CC)CCN1C(CC2=CC=C(OC)C=C2)=NC3=CC=CC=C31',
         'Metodesnitazepyne' : 'COC(C=C1)=CC=C1CC2=NC3=CC=CC=C3N2CCN4CCCC4',
         'Etodesnitazene (Etazene)' : 'CCN(CCN1C(CC2=CC=C(C=C2)OCC)=NC3=CC=CC=C31)CC',
         'Etodesnitazepyne' : 'CCOC(C=C1)=CC=C1CC2=NC3=CC=CC=C3N2CCN4CCCC4',
         'Etodesnitazepipne' : 'CCOC(C=C1)=CC=C1CC2=NC3=CC=CC=C3N2CCN4CCCCC4',
         'Protodesnitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCCC)C=C2)=NC3=CC=CC=C31',
         'Isotodesnitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OC(C)C)C=C2)=NC3=CC=CC=C31',
         'Nitazene' : 'CCN(CCN1C(CC2=CC=CC=C2)=NC3=CC([N+]([O-])=O)=CC=C31)CC',
         'Metonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'meta-Metonitazene' : 'CCN(CC)CCN1C(CC2=CC=CC(OC)=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Metonitazepyne' : 'COC(C=C1)=CC=C1CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCC4',
         'Metonitazepipne' : 'COC(C=C1)=CC=C1CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCCC4',
         'N-Desethylmetonitazene' : 'COC(C=C1)=CC=C1CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCNCC',
         'Metomethazene' : 'CC1=CC=C(N(CCN(CC)CC)C(CC2=CC=C(OC)C=C2)=N3)C3=C1',
         'Dimetonitazene' : 'COC(C(OC)=C1)=CC=C1CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC',
         'α-methyl-metonitazene' : 'COC(C=C1)=CC=C1C(C)C2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC',
         'Metonitazene phenethyl homologue (Ethylene metonitazene)' : 'CCN(CC)CCN1C(CCC2=CC=C(OC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Etonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'O-Desethyl-etonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(O)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'N-Desethyletonitazene (NDE)' : 'CCNCCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Etonitazene 5-amino metabolite' : 'NC1=CC=C(N(CCN(CC)CC)C(CC2=CC=C(OCC)C=C2)=N3)C3=C1',
         'Etomethazene' : 'CC1=CC=C(N(CCN(CC)CC)C(CC2=CC=C(OCC)C=C2)=N3)C3=C1',
         'Etonitazene 5-trifluoromethyl analogue (Etotriflazene)' : 'CCN(CC)CCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC(C(F)(F)F)=CC=C31',
         'Etonitazene_5-cyano_analogue (Etocyanazene)' : 'CCN(CC)CCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC(C#N)=CC=C31',
         'Etonitazene_5-acetyl_analogue (Etoacetazene)' : 'CCN(CC)CCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC(C(C)=O)=CC=C31',
         'Etonitazene 5,6-dichloro analogue (Etodicloazene)' : 'CCN(CC)CCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC(Cl)=C(Cl)C=C31',
         'Etonitazene N,N-dimethyl analogue' : 'CCN(CC)CCN1C(CC2=CC=C(OCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Etonitazepyne' : 'CCOC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCC4)C=C1',
         'Etonitazepipne' : 'CCOC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCCC4)C=C1',
         'Etonitazene morpholine analogue' : 'CCOC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCOCC4)C=C1',
         'Etonitazene 6-nitro isomer' : 'CCOC1=CC=C(CC2=NC3=CC=C([N+]([O-])=O)C=C3N2CCN(CC)CC)C=C1',
         'Protonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Protonitazepyne' : 'CCCOC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCC4)C=C1',
         'Protonitazepipne' : 'CCCOC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCCC4)C=C1',
         'N-Desethylprotonitazene' : 'CCCOC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCNCC)C=C1',
         'Isotonitazene' : 'CC(C)OC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C=C1',
         'Isotonitazepyne' : 'CC(C)OC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCC4)C=C1',
         'Isotonitazepipne' : 'CC(C)OC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN4CCCCC4)C=C1',
         'N-Desethylisotonitazene' : 'CC(C)OC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCNCC)C=C1',
         'Butonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCCCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Isobutonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCC(C)C)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Secbutonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OC(C)CC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Etoetonitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCCOCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Flunitazene' : 'FC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C=C1',
         'Clonitazene' : 'ClC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C=C1',
         'Diclonitazene' : 'ClC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C(Cl)=C1',
         'α-carboxamido-clonitazene' : 'ClC1=CC=C(C(C(N)=O)C2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C=C1',
         'Bronitazene' : 'BrC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C=C1',
         'Nitronitazene' : 'CCN(CC)CCN1C(CC2=CC=C([N+]([O-])=O)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Methylnitazene' : 'CC1=CC=C(CC2=NC3=CC([N+]([O-])=O)=CC=C3N2CCN(CC)CC)C=C1',
         'Ethylnitazene' : 'CCN(CC)CCN1C(CC2=CC=C(CC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Propylnitazene' : 'CCN(CC)CCN1C(CC2=CC=C(CCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         't-Butylnitazene' : 'CCN(CC)CCN1C(CC2=CC=C(C(C)(C)C)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Acetoxynitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OC(C)=O)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Methylthionitazene' : 'CCN(CC)CCN1C(CC2=CC=C(SC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Ethylthionitazene' : 'CCN(CC)CCN1C(CC2=CC=C(SCC)C=C2)=NC3=CC([N+]([O-])=O)=CC=C31',
         'Etodesnitazene phenylthio analogue' : 'CCN(CC)CCN1C(SC2=CC=C(OCC)C=C2)=NC3=CC=CC=C31',
         'Etodesnitazene phenylthio' : 'CCOC1=CC=C(SC2=NC3=CC=CC=C3N2CCN4CCCC4)C=C1',
         'Methylenedioxynitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCO3)C3=C2)=NC4=CC([N+]([O-])=O)=CC=C41',
         'Ethyleneoxynitazene' : 'CCN(CC)CCN1C(CC2=CC=C(OCC3)C3=C2)=NC4=CC([N+]([O-])=O)=CC=C41'
              
         }

         return smiles

    @staticmethod
    def get_smarts():

         smarts = {
         }

         return smarts
