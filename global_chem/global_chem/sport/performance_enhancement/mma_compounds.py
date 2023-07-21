#!/usr/bin/env python3
#
# GlobalChem - MMA COMPOUNDS
#
# -------------------------


class MmaCompounds(object):
    def __init__(self):
        self.name = "mma_compounds"

    @staticmethod
    def get_smiles():
        smiles = {
            "1-Androstenedione": "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CC[C@@H]4[C@@]3(C=CC(=O)C4)C",
            "Dextroamphetamine": "C[C@@H](CC1=CC=CC=C1)N",
            "Anastrozole": "CC(C)(C#N)C1=CC(=CC(=C1)CN2C=NC=N2)C(C)(C)C#N",
            "Arimistane": "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)C(=O)C=C4[C@@]3(CCC=C4)C",
            "Boldenone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)C=C[C@]34C",
            "Cardarine": "CC1=C(C=CC(=C1)SCC2=C(N=C(S2)C3=CC=C(C=C3)C(F)(F)F)C)OCC(=O)O",
            "Clenbuterol": "CC(C)(C)NCC(C1=CC(=C(C(=C1)Cl)N)Cl)O",
            "Clomiphene": "CCN(CC)CCOC1=CC=C(C=C1)/C(=C(\C2=CC=CC=C2)/Cl)/C3=CC=CC=C3",
            "Cocaine": "CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC",
            "Dehydroepiandrosterone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CC=C4[C@@]3(CC[C@@H](C4)O)C",
            "Dehydrochlormethyltestosterone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C)O)CCC4=C(C(=O)C=C[C@]34C)Cl",
            "Furosemide": "C1=COC(=C1)CNC2=CC(=C(C=C2C(=O)O)S(=O)(=O)N)Cl",
            "Hydrochlorothiazide": "C1NC2=CC(=C(C=C2S(=O)(=O)N1)S(=O)(=O)N)Cl",
            "4-amino-6-chloro-1,3-benzenedisulfonamide": "C1=C(C(=CC(=C1Cl)S(=O)(=O)N)S(=O)(=O)N)N",
            "Drostanolone": "C[C@@H]1C[C@]2([C@@H](CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CC[C@@H]4O)C)CC1=O)C",
            "Erythropoietin": "[H]C1=C([H])C(C2=C([H])C([H])=C([H])C([H])=C2C3=NN=NC3)=C([H])C([H])=C1C(N4C(C([H])([H])O)=C(Cl)N=C4C([H])([H])C(C([H])([H])C([H])(C)[H])([H])[H])([H])[H]",
            "GHRP-6": "C[C@@H](C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC3=CC=CC=C3)C(=O)N[C@@H](CCCCN)C(=O)N)NC(=O)[C@@H](CC4=CNC5=CC=CC=C54)NC(=O)[C@H](CC6=CN=CN6)N",
            "GHRP-2": "C[C@H](C(=O)N[C@H](CC1=CC2=CC=CC=C2C=C1)C(=O)N[C@@H](C)C(=O)N[C@@H](CC3=CNC4=CC=CC=C43)C(=O)N[C@H](CC5=CC=CC=C5)C(=O)N[C@@H](CCCCN)C(=O)N)N",
            "Higenamine": "C1CNC(C2=CC(=C(C=C21)O)O)CC3=CC=C(C=C3)O",
            "Ibutamoren": "CC(C)(C(=O)N[C@H](COCC1=CC=CC=C1)C(=O)N2CCC3(CC2)CN(C4=CC=CC=C34)S(=O)(=O)C)N",
            "Ipamorelin": "CC(C)(C(=O)N[C@@H](CC1=CN=CN1)C(=O)N[C@H](CC2=CC3=CC=CC=C3C=C2)C(=O)N[C@H](CC4=CC=CC=C4)C(=O)N[C@@H](CCCCN)C(=O)N)N",
            "Ligandrol": "C1C[C@@H](N(C1)C2=CC(=C(C=C2)C#N)C(F)(F)F)[C@H](C(F)(F)F)O",
            "Meldonium": "C[N+](C)(C)NCCC(=O)[O-]",
            "Methyltestosterone": "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CC[C@]4(C)O)C",
            "Mesterolone": "C[C@H]1CC(=O)C[C@H]2[C@]1([C@H]3CC[C@]4([C@H]([C@@H]3CC2)CC[C@@H]4O)C)C",
            "Methandienone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C)O)CCC4=CC(=O)C=C[C@]34C",
            "Modafinil": "C1=CC=C(C=C1)C(C2=CC=CC=C2)S(=O)CC(=O)N",
            "Nandrolone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@H]34",
            "19-Norandrosterone": "C[C@]12CC[C@@H]3[C@H]4CC[C@H](C[C@@H]4CC[C@H]3[C@@H]1CCC2=O)O",
            "Ostarine": "C[C@](COC1=CC=C(C=C1)C#N)(C(=O)NC2=CC(=C(C=C2)C#N)C(F)(F)F)O",
            "Ozone": "[O-][O+]=O",
            "Stanozolol": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C)O)CC[C@@H]4[C@@]3(CC5=C(C4)NN=C5)C",
            "Tamoxifen": "CC/C(=C(\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3",
            "Exogenous Testosterone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C",
            "Tetrahydrocannabinol": "CC(O1)(C)[C@@]2(C)[C@](C=C(C)CC2)(C)C3=C1C=C(CCCCC)C=C3O",
            "Trenbolone": "C[C@]12C=CC3=C4CCC(=O)C=C4CC[C@H]3[C@@H]1CC[C@@H]2O",
        }

        return smiles
