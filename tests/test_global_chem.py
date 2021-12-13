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
        ,gc.ribose_subpocket_smiles
        ,gc.adenine_subpocket_smiles
        ,gc.hydrophobic_subpocket_smiles
        ,gc.type_1_subpocket_smiles
        ,gc.type_2_subpocket_smiles
        ,gc.exposed_to_solvent_smiles,
        gc.alpha_amino_removed_by_acid_smiles,
        gc.alpha_amino_removed_by_acid_smarts,
        gc.alpha_amino_removed_by_acid_acronym_smiles,
        gc.alpha_amino_removed_by_acid_acronym_smarts,
        gc.alpha_amino_removed_by_base_smiles,
        gc.alpha_amino_removed_by_base_smarts,
        gc.alpha_amino_removed_by_base_acronym_smiles,
        gc.alpha_amino_removed_by_base_acronym_smarts,
        gc.other_alpha_amino_smiles,
        gc.other_alpha_amino_smarts,
        gc.other_alpha_amino_acronym_smiles,
        gc.other_alpha_amino_acronym_smarts,
        gc.lys_orn_dap_dab_removed_by_acid_smiles,
        gc.lys_orn_dap_dab_removed_by_acid_smarts,
        gc.lys_orn_dap_dab_removed_by_acid_acronym_smiles,
        gc.lys_orn_dap_dab_removed_by_acid_acronym_smarts,
        gc.lys_orn_dap_dab_removed_by_base_smiles,
        gc.lys_orn_dap_dab_removed_by_base_smarts,
        gc.lys_orn_dap_dab_removed_by_base_acronym_smiles,
        gc.lys_orn_dap_dab_removed_by_base_acronym_smarts,
        gc.other_lys_orn_dap_dab_smiles,
        gc.other_lys_orn_dap_dab_smarts,
        gc.other_lys_orn_dap_dab_acronym_smiles,
        gc.other_lys_orn_dap_dab_acronym_smarts,
        gc.alpha_carboxylic_acid_removed_by_acid_smiles,
        gc.alpha_carboxylic_acid_removed_by_acid_smarts,
        gc.alpha_carboxylic_acid_removed_by_acid_acronym_smiles,
        gc.alpha_carboxylic_acid_removed_by_acid_acronym_smarts,
        gc.alpha_carboxylic_acid_removed_by_base_smiles,
        gc.alpha_carboxylic_acid_removed_by_base_smarts,
        gc.alpha_carboxylic_acid_removed_by_base_acronym_smiles,
        gc.alpha_carboxylic_acid_removed_by_base_acronym_smarts,
        gc.other_alpha_carboxylic_acid_protecting_group_smiles,
        gc.other_alpha_carboxylic_acid_protecting_group_smarts,
        gc.other_alpha_carboxylic_acid_protecting_group_acronym_smiles,
        gc.other_alpha_carboxylic_acid_protecting_group_acronym_smarts,
        gc.asp_glu_removed_by_acid_smiles,
        gc.asp_glu_removed_by_acid_smarts,
        gc.asp_glu_removed_by_acid_acronym_smiles,
        gc.asp_glu_removed_by_acid_acronym_smarts,
        gc.asp_glu_removed_by_base_smiles,
        gc.asp_glu_removed_by_base_smarts,
        gc.asp_glu_removed_by_base_acronym_smiles,
        gc.asp_glu_removed_by_base_acronym_smarts,
        gc.other_asp_glu_smiles,
        gc.other_asp_glu_smarts,
        gc.other_asp_glu_acronym_smiles,
        gc.other_asp_glu_acronym_smarts,
        gc.amide_backbone_protecting_group_removed_by_acid_smiles,
        gc.amide_backbone_protecting_group_removed_by_acid_smarts,
        gc.amide_backbone_protecting_group_removed_by_acid_acronym_smiles,
        gc.amide_backbone_protecting_group_removed_by_acid_acronym_smarts,
        gc.other_amide_backbone_protecting_group_smiles,
        gc.other_amide_backbone_protecting_group_smarts,
        gc.asn_gln_removed_by_acid_smiles,
        gc.asn_gln_removed_by_acid_smarts,
        gc.asn_gln_removed_by_acid_acronym_smiles,
        gc.asn_gln_removed_by_acid_acronym_smarts,
        gc.arg_removed_by_acid_smiles,
        gc.arg_removed_by_acid_smarts,
        gc.arg_removed_by_acid_acronym_smiles,
        gc.arg_removed_by_acid_acronym_smarts,
        gc.arg_removed_by_base_smiles,
        gc.arg_removed_by_base_smarts,
        gc.arg_removed_by_base_acronym_smiles,
        gc.arg_removed_by_base_acronym_smarts,
        gc.other_arg_smiles,
        gc.other_arg_smarts,
        gc.other_arg_acronym_smiles,
        gc.other_arg_acronym_smarts,
        gc.cys_removed_by_acid_smiles,
        gc.cys_removed_by_acid_smarts,
        gc.cys_removed_by_acid_acronym_smiles,
        gc.cys_removed_by_acid_acronym_smarts,
        gc.cys_removed_by_base_smiles,
        gc.cys_removed_by_base_smarts,
        gc.cys_removed_by_base_acronym_smiles,
        gc.cys_removed_by_base_acronym_smarts,
        gc.other_cys_smiles,
        gc.other_cys_smarts,
        gc.other_cys_acronym_smiles,
        gc.other_cys_acronym_smarts,
        gc.his_removed_by_acid_smiles,
        gc.his_removed_by_acid_smarts,
        gc.his_removed_by_acid_acronym_smiles,
        gc.his_removed_by_acid_acronym_smarts,
        gc.his_removed_by_base_smiles,
        gc.his_removed_by_base_smarts,
        gc.his_removed_by_base_acronym_smiles,
        gc.his_removed_by_base_acronym_smarts,
        gc.other_his_smiles,
        gc.other_his_smarts,
        gc.other_his_acronym_smiles,
        gc.other_his_acronym_smarts,
        gc.ser_thr_hyp_removed_by_acid_smiles,
        gc.ser_thr_hyp_removed_by_acid_smarts,
        gc.ser_thr_hyp_removed_by_acid_acronym_smiles,
        gc.ser_thr_hyp_acronym_removed_by_acid_smarts,
        gc.other_ser_thr_hyp_smiles,
        gc.other_ser_thr_hyp_smarts,
        gc.other_ser_thr_hyp_acronym_smiles,
        gc.other_ser_thr_hyp_acronym_smarts,
        gc.tyr_removed_by_acid_smiles,
        gc.tyr_removed_by_acid_smarts,
        gc.tyr_removed_by_acid_acronym_smiles,
        gc.tyr_removed_by_acid_acronym_smarts,
        gc.other_tyr_protecting_group_smiles,
        gc.other_tyr_protecting_group_smarts,
        gc.other_tyr_protecting_group_acronym_smiles,
        gc.other_tyr_protecting_group_acronym_smarts,
        gc.trp_removed_by_acid_smiles,
        gc.trp_removed_by_acid_smarts,
        gc.trp_removed_by_acid_acronym_smiles,
        gc.trp_removed_by_acid_acronym_smarts,
        gc.polyfluoroalkyl_smiles,
        gc.montmorillonite_chemical_adsorption_smiles
    ]

    total_molecules = []
    passing_molecules = []

    for i in range(0, len(compounds)):

        molecules = compounds[i]
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

