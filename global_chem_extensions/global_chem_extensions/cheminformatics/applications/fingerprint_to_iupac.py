#!/usr/bin/env python3
#
# GlobalChemExtensions - Classify Chemical Fingerprints
#
# ------------------------------------------------------

class ClassifyChemicalFingerprints(object):

    '''

    Return a classification layer to RDKF Fingerprints for Daylight Substructure

    '''

    def __init__(self, bit_info = {}, smiles = None):

        self.bit_info = bit_info
        self.smiles = smiles
        self.smiles_fragments = []

        self.fingerprint = None

    def generate_fingerprint(self):

        fingerprint = Chem.RDKFingerprint(
            molecule,
            bitInfo=self.bit_info,
            maxPath=7,
            minPath=6
        )

        self.fingerprint = fingerprint

    def generate_smiles_map(self):


        '''

        Convert the fingerprint to IUPAC name based on the GlobalChem

        '''

        for key, bond_indices in self.bit_info.items():

            atoms = set()

            for bond_index in bond_indices[0]:

                atoms.add(molecule.GetBondWithIdx(bond_index).GetBeginAtomIdx())
                atoms.add(molecule.GetBondWithIdx(bond_index).GetEndAtomIdx())

            self.smiles_fragments.append(
                Chem.MolFragmentToSmiles(
                    molecule,
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
                print (all_names[smiles_index])

            except Exception as e:
                continue

        return iupac_map

#
# if __name__ == '__main__':
#
#     from rdkit import Chem
#
#     molecule = Chem.MolFromSmiles('CCOc1ncccc1')
#
#     fingerprint_classifier = ClassifyChemicalFingerprints(
#         smiles = 'CCOc1ncccc1'
#     )
#
#     fingerprint_classifier.generate_fingerprint()
#     smiles_map = fingerprint_classifier.generate_smiles_map()
#     iupac_map = fingerprint_classifier.generate_iupac_map()
#
#     print (iupac_map)