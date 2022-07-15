#!/usr/bin/env python3
#
# GlobalChemExtensions - Classify Chemical Fingerprints
#
# ------------------------------------------------------

# Imports
# -------

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

from global_chem import GlobalChem

functional_groups = {
    'tabun': 'CCOP(=O)(C#N)N(C)C',
    'sarin': 'CC(C)OP(=O)(C)F',
    'soman': 'CC(C(C)(C)C)OP(=O)(C)F',
    'cyclosarin': 'CP(=O)(OC1CCCCC1)F',
    'vx': 'CCOP(C)(=O)SCCN(C(C)C)C(C)C',
    'russian vx': 'CCN(CC)CCSP(=O)(C)OCC(C)C',
    'mirzayanov-a230': 'CCN(CC)C(C)=N[P](C)(F)=O',
    'mirzayanov-a232': 'CCN(CC)C(\C)=N\P(F)(=O)OC',
    'mirzayanov-a234': r'CCOP(F)(=O)\N=C(/C)\N(CC)CC',
    'hoenig-a230': r'Cl/C(F)=N/OP(F)(OCCCl)=O',
    'hoenig-a232': r'Cl/C(F)=N/OP(F)(OC(C)CCl)=O',
    'hoenig-a234': r'Cl/C(F)=N/OP(F)(OC(C)C(C)Cl)=O',
    'novichok-5': 'FP1OC(C)CO1',
    'novichok-7': 'FP1OC(C)C(C)O1',
}

class DecoderEngine(object):

    def __init__(self, fingerprint = None):

        self.fingerprint = fingerprint
        self.gc = GlobalChem()

    @staticmethod
    def generate_morgan_fingerprint(smiles):

        '''

        Generate the Morgan Fingerprint

        Defaults:
            Bit_representation: needs to be 512 to capture hollistic information about the structure and not as acute.
            Can be modified later.

        '''

        molecule = Chem.MolFromSmiles(smiles)

        bit_string = AllChem.GetMorganFingerprintAsBitVect(
            molecule,
            2,
            nBits=512
        ).ToBitString()

        return bit_string

    def classify_fingerprint(self, fingerprint, node='cengage'):

        scores = []
        identified_functional_groups = []
        all_smiles = []

        fingerprint = DataStructs.cDataStructs.CreateFromBitString(fingerprint)

        for smiles, reference_fingerprint in bit_strings.items():

            reference_fingerprint = DataStructs.cDataStructs.CreateFromBitString(reference_fingerprint)
            score = DataStructs.FingerprintSimilarity(fingerprint, reference_fingerprint)

            all_smiles.append(smiles)
            scores.append(score)

        for i, score in enumerate(scores):

            if score > 0.90:
                identified_functional_groups.append(all_smiles[i])

        print (identified_functional_groups)

    def classify_smiles(self, node='all'):

        '''

        Classify the SMILES of a molecule by decomposing into it's functional groups
        using the BRICS module and then generating fingerprints.

        '''



if __name__ == '__main__':

    from rdkit import Chem
    from rdkit.Chem import BRICS

    # import re
    #
    # fragments = list(BRICS.BRICSDecompose(Chem.MolFromSmiles('CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3')))
    #
    # pattern_match_1 = re.compile(r'(?!\])\(\[\d+\*]\).*?', flags=re.MULTILINE)
    # pattern_match_2 = re.compile(r'(?!\])\[\d+\*].*?', flags=re.MULTILINE)
    #
    # cleaned_fragments = []
    #
    # for fragment in fragments:
    #
    #     match_1 = pattern_match_1.findall(fragment)
    #     match_2 = pattern_match_2.findall(fragment)
    #
    #     if match_1:
    #         for i in match_1:
    #             fragment = fragment.replace(i, '')
    #
    #     if match_2:
    #         for i in match_2:
    #             fragment = fragment.replace(i, '')
    #
    #     cleaned_fragments.append(fragment)
    #
    decoder_engine = DecoderEngine()
    #
    # for cleaned_fragment in cleaned_fragments:
    #
    #     morgan_fingerprint  = decoder_engine.generate_morgan_fingerprint(
    #         cleaned_fragment
    #     )
    #
    #     decoder_engine.classify_fingerprint(
    #         morgan_fingerprint
    #     )

    for k, v in functional_groups.items():
        try:
            print (f"'{k}': '{decoder_engine.generate_morgan_fingerprint(v)}',")
        except:
            continue

