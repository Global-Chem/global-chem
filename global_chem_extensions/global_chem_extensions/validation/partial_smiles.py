#!/usr/bin/env python3
#
# GlobalChemExtensions - Partial SMILES
#
# --------------------------------------

# Imports
# -------

import partialsmiles as ps
from pysmiles import read_smiles
from molvs import validate_smiles

# RDKit Imports
# -------------

from rdkit import Chem

class PartialSmilesValidation(object):

    __version__ = '0.0.1'


    def __init__(self,
                 partial=False
                 ):

        self.partial = partial

    def validate(self, smiles_list, partial_smiles=True, rdkit=True, pysmiles=True, molvs=True):

        '''

        Validate the SMILES

        Arguments:
            smiles_list (List): List of SMILES that the user would like to put
            partial_smiles (Bool): Partial SMILES
            rdkit (Bool): whether to include rdkit
            pysmiles (Bool): whether to include pysmiles
            molvs (Bool): whether to include molvs

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
                if pysmiles:
                    mol = read_smiles(smiles)
                if molvs:
                    errors = validate_smiles(smiles)
                    print (errors)
                successes.append(smiles)
            except:
                failures.append(smiles)

        return successes, failures