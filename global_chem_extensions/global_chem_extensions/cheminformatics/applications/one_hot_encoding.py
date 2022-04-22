#!/usr/bin/env python3
#
# GlobalChemExtensions - One-Hot Encoder

# --------------------------------------

# Imports
# -------

import numpy as np


class SmilesOneHotEncoder(object):

    __version__ = '0.0.1'

    __SMILES_MAPPING__ = [
        ' ',
        '#', '%', '(', ')', '+', '-', '.', '/',
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        '=', '@',
        'A', 'B', 'C', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P',
        'R', 'S', 'T', 'V', 'X', 'Z',
        '[', '\\', ']',
        'a', 'b', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'p', 'r', 's',
        't', 'u',
        '&', ':', '*'
    ]

    def __init__(
            self,
            smiles_list,
            max_length = 120
        ):

        self.smiles_list = smiles_list
        self.max_length = max_length

        self.smiles_to_index = None
        self.index_to_smiles = None

        self._setup_encoders()

    def _setup_encoders(self):

        '''

        Setup the conversion encoders between the SMILES language and the Indexes

        '''

        self.smiles_to_index = dict( (c,i) for i,c in enumerate( self.__SMILES_MAPPING__ ) )
        self.index_to_smiles = dict( (i,c) for i,c in enumerate( self.__SMILES_MAPPING__ ) )

    def encode(self):

        '''

        Encode the List of SMILES

        Returns:

            encoded_smiles_list (List): One Hot Encoded of the SMILES List.

        '''

        encoded_smiles_list = []

        for smiles in self.smiles_list:

            encoded_smiles = np.zeros(
                (
                    self.max_length,
                    len( self.__SMILES_MAPPING__ )
                )
            )

            for i, c in enumerate( smiles ):

                encoded_smiles[i, self.smiles_to_index[c] ] = 1

            encoded_smiles_list.append(encoded_smiles)

        return encoded_smiles_list

    def decode(self):


        '''

        Decode the SMILES from the index

        Returns:

            decoded_smiles (List): List of the decoded SMILES from the index mapping

        '''

        decoded_smiles = []

        for encoded_smiles in self.smiles_list:

            smiles = ''

            encoded_smiles = encoded_smiles.argmax( axis=-1 )

            for i in encoded_smiles:
                smiles += self.index_to_smiles[ i ]

            decoded_smiles.append(smiles)

        return decoded_smiles