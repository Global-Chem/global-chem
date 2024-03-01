#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Lanthipeptides

# -------------------------------------------------

class GlobalChemLanthiPeptide(object):

    __AMINO_ACID_SEQUENCE__ = {
        "A": "C",
        "R": "CCCCNC(N)=N",
        "N": "CC(N)=O",
        "D": "CC(O)=O",
        "B": "CC(O)=O",
        "C": "CS",
        "E": "CCC(O)=O",
        "Q": "CCC(N)=O",
        "Z": "CCC(N)=O",
        "G": "[H]",
        "H": "CC1=CNC=N1",
        "I": "C(CC)([H])C",
        "L": "CC(C)C",
        "K": "CCCCN",
        "M": "CCSC",
        "F": "CC1=CC=CC=C1",
        "P": "C2CCCN2",
        "S": "CO",
        "T": "C(C)([H])O",
        "W": "CCC1=CNC2=C1C=CC=C2",
        "Y": "CC1=CC=C(O)C=C1",
        "V": "C(C)C"
    }

    def __init__(self, name = '', length = 1):

        self.name = name
        self.length = length



if __name__ == '__main__':

    nisin = GlobalChemLanthiPeptide(
        name='nisin'
    )

    mutacin = GlobalChemLanthiPeptide(
        name='mutacin'
    )

    print (nisin.name)
    print (mutacin.name)