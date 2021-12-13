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
        gc.get_amino_acids(),
        gc.get_rings_in_drugs(),
        gc.get_essential_vitamins(),
        gc.get_common_organic_solvents(),
        gc.get_polyfluoroalkyl_substances(),
        gc.get_common_privileged_scaffolds(),
        gc.get_open_smiles_functional_groups(),
        gc.get_common_polymer_repeating_units(),
        gc.get_braf_kinase_inhibitors_for_cancer(),
        gc.get_common_heterocyclic_rings_phase_2(),
        gc.get_commonly_used_r_group_replacements(),
        gc.get_common_warhead_covalent_inhibitors(),
        gc.get_common_amino_acid_protecting_groups(),
        gc.get_iupac_blue_book_common_functional_groups(),
        gc.get_common_electrophilic_warheads_for_kinases(),
        gc.get_privileged_scaffolds_for_kinase_inhibitors(),
        gc.get_chemical_adsorption_on_montmorillonite_clays(),
    ]

    total_molecules = []
    passing_molecules = []

    for i in range(0, len(compounds)):

        molecule_set = list(compounds[i])

        for molecules in molecule_set:
            smiles_list = list(list(molecules.values()))

            for i in range(0, len(smiles_list)):

                total_molecules.append(smiles_list[i])

                try:
                    mol = Chem.MolFromSmiles(smiles_list[i])
                    if mol:
                        passing_molecules.append(smiles_list[i])
                except Exception as e:
                    pass

    print ('Total Molecules: %s ' % (len(total_molecules)))
    print ('Total RDKit Passing Molecules: %s ' % len(passing_molecules))
