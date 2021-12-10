# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

# imports
# -------

from rdkit import Chem
from global_chem import GlobalChem

def test_rdkit_passing():

    '''

    Test the global chem class through rdkit parser

    '''

    gc = GlobalChem()

    compounds = [
        gc.amino_acid_side_smiles
        ,gc.common_organic_solvents_smiles
        ,gc.common_polymer_repeating_units_smiles
        ,gc.common_warheads_smiles
        ,gc.iupac_blue_book_radical_smiles
        ,gc.iupac_blue_book_rings_smiles
        ,gc.kinase_electrophilic_warheads_smiles
        ,gc.kinase_privileged_scaffolds_smiles
        ,gc.open_smiles_functional_groups_smiles
        ,gc.phase_2_hetereocyclic_rings_smiles
        ,gc.privileged_scaffolds_smiles
        ,gc.r_groups_replacements_smiles
        ,gc.rings_in_drugs_smiles
        ,gc.vitamins_smiles
    ]

    total_molecules = []
    passing_rdkit_molecules = []

    for i in range(0, len(compounds)):

        molecules = compounds[i]

        smiles_list = list(list(molecules.values()))

        for i in range(0, len(smiles_list)):

            total_molecules.append(smiles_list[i])

            try:
                mol = Chem.MolFromSmiles(smiles_list[i])
                passing_rdkit_molecules.append(smiles_list[i])

            except:
                pass

    print ('Total Molecules: %s ' % len(total_molecules))
    print ('Total Passing RDKit: %s ' % len(passing_rdkit_molecules))

