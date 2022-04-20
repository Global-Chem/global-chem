#!/usr/bin/env python3
#
# GlobalChemExtensions - Partial SMILES
#
# --------------------------------------

# Imports
# -------

import urllib.request
from urllib.parse import quote

# Validation Imports
# ------------------

import deepsmiles
import selfies as sf
import partialsmiles as ps

from rdkit import Chem
from pysmiles import read_smiles
from molvs import validate_smiles

class PartialSmilesValidation(object):

    __version__ = '0.0.1'


    def __init__(self,
                 partial=False
                 ):

        self.partial = partial
        self.deep_smiles_converter = deepsmiles.Converter(rings=True, branches=True)

    def validate(
            self,
            smiles_list,
            partial_smiles = True,
            rdkit = False,
            pysmiles = False,
            molvs = False,
            deepsmiles = False,
            selfies = False,
    ):

        '''

        Validate the SMILES

        Arguments:
            smiles_list (List): List of SMILES that the user would like to put
            partial_smiles (Bool): Partial SMILES
            rdkit (Bool): whether to include rdkit
            pysmiles (Bool): whether to include pysmiles
            molvs (Bool): whether to include molvs
            deepsmiles (Bool): deepSMILES validation for machine learning
            selfies (Bool): SELFIES validation for machine learning

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
                    for error in errors:
                        if 'ERROR' in error:
                            raise Exception
                if deepsmiles:

                    encoded = self.deep_smiles_converter.encode(smiles)
                    try:
                        encoded = self.deep_smiles_converter.encode(smiles)
                        decoded = self.deep_smiles_converter.decode(encoded)
                    except deepsmiles.EncodeError as e:
                        continue
                    except deepsmiles.DecodeError as e:
                        continue

                if selfies:
                    try:
                        encoded = sf.encoder(smiles)
                        decoded = sf.decoder(encoded)
                    except sf.EncoderError:
                        pass
                    except sf.DecoderError:
                        pass
                successes.append(smiles)

            except Exception as e:
                failures.append(smiles)

        return successes, failures