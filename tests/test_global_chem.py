# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

# imports
# -------
from indigo import *
from rdkit import Chem
from global_chem import GlobalChem

def test_rdkit_passing():

    '''

    Test the global chem class through rdkit parser

    '''

    total_smiles = GlobalChem().get_all_smiles()

    passing_molecules = []

    for i in range(0, len(total_smiles)):

        try:
            mol = Chem.MolFromSmiles(total_smiles[i])
            if mol:
                passing_molecules.append(total_smiles[i])
        except Exception as e:
            pass

    print ('Total Molecules: %s ' % (len(total_smiles)))
    print ('Total RDKit Passing Molecules: %s ' % len(passing_molecules))
    print ('Failed Compounds: %s' % list(set(total_smiles) - set(passing_molecules)))

def test_indigo_passing():

    '''

    Test the GlobalChem Passing through the indigo parser.

    '''

    total_smiles = GlobalChem().get_all_smiles()

    indigo = Indigo()

    success_compounds = []
    failed_compounds = []

    for smiles in total_smiles:

        try:
            molecule = indigo.loadMolecule(smiles)
            success_compounds.append(smiles)
        except IndigoException as e:
            failed_compounds.append(smiles)

    print ('Total Molecules: %s ' % (len(total_smiles)))
    print ('Total Indigo Passing Molecules: %s ' % len(success_compounds))
    print ('Failed Compounds: %s' % len(failed_compounds))
    print ('Failed Compounds: %s' % failed_compounds)