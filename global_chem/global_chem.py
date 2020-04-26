#!/usr/bin/env python3
#
# GlobalChem - Content Variable Store
#
# -----------------------------------


class GlobalChem(object):

    __version__ = "0.1.0"
    __allow_update__ = True

    """

    GlobalChem will be the master class of all variables, as the content store grows we can use this as the parent class.

    """


    def _get_amino_acids(self):

        amino_acid_side_chains = {
            "alanine": "C",  "arginine": "CCCCNC(N)=N", "asparagine": "CCC(N)=O", "aspartic acid": "CC(O)=O",
            "cysteine": "CS", "glutamic acid": "CCC(O)=O", "glutamine": "CCC(N)=O", "glycine": "[H]",
            "histidine": "CC1=CNC=N1", "isoleucine": "C(CC)([H])C", "leucine": "CC(C)C", "lysine": "CCCCN",
            "methionine": "CCSC", "phenylalanine": "CC1=CC=CC=C1", "proline": "C2CCCN2", "serine": "CO",
            "threonine": "C(C)([H])O", "tryptophan": "CCC1=CNC2=C1C=CC=C2", "tyrosine": "CC1=CC=C(O)C=C1",
            "valine": "C(C)C"
        }

        return amino_acid_side_chains

    def _get_functional_groups_smiles(self):

        functional_groups_smies = {
            "bromine": "Br", "chlorine": "Cl", "fluorine": "F", "acyl bromide": "C(=O)Br",
            "acyl chloride": "C(=O)Cl", "acyl fluoride": "C(=O)F", "acyl iodide": "C(=O)I",
            "primary alcohol": "O", "ketone": "C(=O)OC", "carboxylic acid": "C(=O)O",
            "acetic anydride": "CC(=O)OC(=O)C", "primary amine": "N",
            "secondary amine": "NC", "enamine": "N", "amide": "C(=O)N",
            "nitro": "[N+](=O)[O-]", "sulfoxide": "S(=O)(=O)", "ether": "COC", "azide": "C([N-][N+]#N)"
        }

        return functional_groups_smies

    def _get_functional_groups_smarts(self):

        functional_groups_smarts = {
            "bromine": "[Br]", "chlorine": "[Cl]", "fluorine": "[F]", "acyl bromide": "[CX3](=[OX1])[Br]",
            "acyl chloride": "[CX3](=[OX1])[Cl]", "acyl fluoride": "[CX3](=[OX1])[F]",
            "acyl iodide": "[CX3](=[OX1])[I]", "primary alcohol": "[OX2H]", "ketone": "[CX3]=[OX1]",
            "carboxylic acid": "[CX3](=O)[OX1H0-,OX2H1]", "acetic anydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
            "primary amine": "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]", "secondary amine": "[NX3;H2,H1;!$(NC=O)]",
            "enamine": "[NX3][CX3]=[CX3]", "amide": "[NX3][CX3](=[OX1])[#6]", "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
            "sulfoxide": "[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]",
            "ether": "[OD2]([#6])[#6]", "azide": "[$(-[NX2-]-[NX2+]#[NX1]),$(-[NX2]=[NX2+]=[NX1-])]"
        }

        return functional_groups_smarts

    #------------------------- Property Declaration for GlobalChem ---------------------------#

    amino_acid_side_chains = property(_get_amino_acids)
    functional_groups_smiles = property(_get_functional_groups_smiles)
    functional_groups_smarts = property(_get_functional_groups_smarts)