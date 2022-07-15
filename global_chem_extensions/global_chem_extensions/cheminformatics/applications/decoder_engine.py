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

class ClassifyChemicalFingerprints(object):

    '''

    Return a classification layer to RDKF Fingerprints for Daylight Substructure

    '''

    def __init__(self, bit_info = {}, smiles = None):

        self.bit_info = bit_info
        self.smiles = smiles
        self.molecule = Chem.MolFromSmiles(self.smiles)
        self.smiles_fragments = []

        # Generate fingerprints

        self.fingerprint = self.generate_morgan_fingerprint(self.molecule)

    def generate_fragments(self, molecule):


        '''

        Generate the Fragments of the molecule using the BRICS Module

        '''



    def generate_fingerprint(self):

        fingerprint = Chem.RDKFingerprint(
            self.molecule,
            bitInfo=self.bit_info,
            maxPath=10,
            minPath=3
        )

        self.fingerprint = fingerprint

        return fingerprint

    @staticmethod
    def generate_morgan_fingerprint(molecule):

        '''

        Generate the Morgan Fingerprint

        Defaults:
            Bit_representation: needs to be 512 to capture hollistic information about the structure and not as acute.
            Can be modified later.

        '''

        dummy_array = np.zeros((0,))

        fingerprint = AllChem.GetMorganFingerprintAsBitVect(
            molecule,
            1,
            nBits=512
        )
        DataStructs.ConvertToNumpyArray(fingerprint, dummy_array)

        return dummy_array

    def classify_globalchem_fingerprint(self):

        '''

        Compare fingerprint arrays against each other for most probably match

        '''

        from global_chem import GlobalChem


        gc = GlobalChem()

        all_names = gc.get_all_names()
        all_smiles = gc.get_all_smiles()

        all_molecules = [
            Chem.MolFromSmiles(smiles) for smiles in all_smiles
        ]

        all_molecules = []

    def generate_smiles_map(self):


        '''

        Convert the fingerprint to IUPAC name based on the GlobalChem

        '''

        for key, bond_indices in self.bit_info.items():

            atoms = set()

            for bond_index in bond_indices[0]:

                atoms.add(self.molecule.GetBondWithIdx(bond_index).GetBeginAtomIdx())
                atoms.add(self.molecule.GetBondWithIdx(bond_index).GetEndAtomIdx())

            self.smiles_fragments.append(
                Chem.MolFragmentToSmiles(
                    self.molecule,
                    atomsToUse=list(atoms),
                    bondsToUse=self.bit_info[key][0])
            )

        self.smiles_fragments = list(set(self.smiles_fragments))

        return self.smiles_fragments

    def generate_iupac_map(self):


        '''

        Generate the IUPAC map from the SMILES string

        '''


        from global_chem import GlobalChem

        iupac_map = []

        gc = GlobalChem()

        all_names = gc.get_all_names()
        all_smiles = gc.get_all_smiles()

        lower_all_smiles = [i.lower() for i in all_smiles]

        for fragment in self.smiles_fragments:

            try:

                smiles_index = lower_all_smiles.index(fragment)
                print (all_smiles[smiles_index])

            except Exception as e:
                continue

        return iupac_map


if __name__ == '__main__':

    from rdkit import Chem

    fingerprint_classifier = ClassifyChemicalFingerprints(
        smiles = 'C1=CC=CC=N1'
    )

    fingerprint_classifier = ClassifyChemicalFingerprints(
        smiles = 'c1ccncc1'
    )

    # fingerprint_classifier.classify_globalchem_fingerprint()

    # print(fingerprint_classifier.generate_morgan_fingerprint())

    fingerprint = fingerprint_classifier.generate_fingerprint()
    smiles_map = fingerprint_classifier.generate_smiles_map()
    # print (fingerprint)

    # iupac_map = fingerprint_classifier.generate_iupac_map()
    #
    # print (iupac_map)