#!/usr/bin/env python3
#
# GlobalChemExtensions - Partial SMILES
#
# --------------------------------------

# Imports
# -------

import partialsmiles as ps

# RDKit Imports
# -------------

from rdkit import Chem

class PartialSmilesValidation(object):

    __version__ = '0.0.1'


    def __init__(self,
                 partial=False
                 ):

        self.partial = partial

    def validate(self, smiles_list, partial_smiles=True, rdkit=True):

        '''

        Validate the SMILES

        Arguments:
            smiles_list (List): List of SMILES that the user would like to put

        Returns:
            successes (List): List of SMILES that worked
            failures (List): List of failed smiles that didn't work.

        '''


        successes = []
        failures = []

        for smiles in smiles_list:
            try:
                if partial_smiles:
                    mol = ps.ParseSmiles(smiles, partial=self.partial)
                if rdkit:
                    mol = Chem.MolFromSmiles(smiles)
                successes.append(smiles)
            except:
                failures.append(smiles)



        return successes, failures