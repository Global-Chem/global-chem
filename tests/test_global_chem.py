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

    protecting = []

    for i in list(gc.get_common_amino_acid_protecting_groups()):

        for j in i.values():

            if '[#6' not in j:
                protecting.append(j)

    total_smiles = []

    for i in range(0, len(compounds)):

        molecule_set = list(compounds[i])

        for molecules in molecule_set:
            names = list(list(molecules.keys()))
            smiles_list = list(list(molecules.values()))

            for i in range(0, len(smiles_list)):

                if '[#6]' not in smiles_list[i]:

                    total_smiles.append(smiles_list[i])

    for i in protecting:
        total_smiles.append(i)

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

    protecting = []
    names = []

    for i in list(gc.get_common_amino_acid_protecting_groups()):

        names_1 = list(i.keys())
        smiles = list(i.values())

        for j in range(0, len(smiles)):

            if '[#6' not in smiles[j]:
                protecting.append(smiles[j])
                names.append(names_1[j])

    total_smiles = []

    for i in range(0, len(compounds)):

        molecule_set = list(compounds[i])

        for molecules in molecule_set:
            names = list(list(molecules.keys()))
            smiles_list = list(list(molecules.values()))

            for i in range(0, len(smiles_list)):

                if '[#6]' not in smiles_list[i] and 'X' not in smiles_list[i]:

                    total_smiles.append(smiles_list[i])


    for i in protecting:
        if '[#6]' not in i and 'X' not in i:
            total_smiles.append(i)

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

