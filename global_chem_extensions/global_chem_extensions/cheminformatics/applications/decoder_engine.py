#!/usr/bin/env python3
#
# GlobalChemExtensions - Classify Chemical Fingerprints
#
# ------------------------------------------------------

# Imports
# -------

import re
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import BRICS

from global_chem import GlobalChem

class DecoderEngine(object):

    def __init__(self, fingerprint = None):

        self.fingerprint = fingerprint
        self.gc = GlobalChem()
        self.gc.build_global_chem_network()

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

    def classify_fingerprint(self, fingerprint, node='organic_and_inorganic_bronsted_acids'):

        scores = []
        identified_functional_groups = []
        all_smiles = []

        if node == 'all':
            bit_strings = self.gc.get_all_bits()
        else:
            bit_strings = self.gc.get_node_bits(node)

        fingerprint = DataStructs.cDataStructs.CreateFromBitString(fingerprint)

        for smiles, reference_fingerprint in bit_strings.items():

            reference_fingerprint = DataStructs.cDataStructs.CreateFromBitString(reference_fingerprint)
            score = DataStructs.FingerprintSimilarity(fingerprint, reference_fingerprint)

            all_smiles.append(smiles)
            scores.append(score)

        for i, score in enumerate(scores):

            if score > 0.90:
                identified_functional_groups.append(all_smiles[i])

        return identified_functional_groups

    def classify_smiles_using_bits(self, smiles, node='organic_and_inorganic_bronsted_acids'):

        '''

        Classify the SMILES of a molecule by decomposing into it's functional groups
        using the BRICS module and then generating fingerprints.

        '''

        fragments = list(BRICS.BRICSDecompose(Chem.MolFromSmiles(smiles)))

        pattern_match_1 = re.compile(r'(?!\])\(\[\d+\*]\).*?', flags=re.MULTILINE)
        pattern_match_2 = re.compile(r'(?!\])\[\d+\*].*?', flags=re.MULTILINE)

        cleaned_fragments = []

        for fragment in fragments:

            match_1 = pattern_match_1.findall(fragment)
            match_2 = pattern_match_2.findall(fragment)

            if match_1:
                for i in match_1:
                    fragment = fragment.replace(i, '')

            if match_2:
                for i in match_2:
                    fragment = fragment.replace(i, '')

            cleaned_fragments.append(fragment)

        decoder_engine = DecoderEngine()

        functional_groups = []

        for cleaned_fragment in cleaned_fragments:

            morgan_fingerprint  = decoder_engine.generate_morgan_fingerprint(
                cleaned_fragment
            )

            classified_fingerprint = decoder_engine.classify_fingerprint(
                morgan_fingerprint,
                node=node
            )

            functional_groups.append(classified_fingerprint)

        functional_groups = sum(functional_groups, [])

        return functional_groups

