#!/usr/bin/env python3
#
# GlobalChem - Content Variable Store
#
# -----------------------------------

from .environment.emerging_perfluoroalkyls import EmergingPerFluoroAlkyls

from .materials.clay.montmorillonite_adsorption import MontmorilloniteAdsorption
from .materials.polymers.common_monomer_repeating_units import CommonMonomerRepeatingUnits

from .medicinal_chemistry.warheads.electrophillic_warheads_for_kinases import ElectrophilicWarheadsForKinases
from .medicinal_chemistry.warheads.common_warheads_covalent_inhibitors import CommonWarheadsCovalentInhibitors

from .medicinal_chemistry.rings.rings_in_drugs import RingsInDrugs
from .medicinal_chemistry.rings.iupac_blue_book_rings import IUPACBlueBookRings
from .medicinal_chemistry.rings.phase_2_hetereocyclic_rings import Phase2HetereoCyclicRings

from .medicinal_chemistry.scaffolds.privileged_scaffolds import PrivilegedScaffolds
from .medicinal_chemistry.scaffolds.iupac_blue_book_substituents import IUPACBlueBook
from .medicinal_chemistry.scaffolds.common_r_group_replacements import CommonRGroupReplacements

from .proteins.kinases.braf.common_groups_by_pocket import BRAFInhibitorsByPocket
from .proteins.kinases.scaffolds.privileged_scaffolds import PrivilegedKinaseInhibitorScaffolds

from global_chem.organic_synthesis.solvents.common_organic_solvents import CommonOrganicSolvents
from global_chem.organic_synthesis.protecting_groups.amino_acid_protecting_groups import AminoAcidProtectingGroups

from .miscellaneous.vitamins import Vitamins
from .miscellaneous.open_smiles import OpenSmiles
from .miscellaneous.amino_acids import AminoAcids
from .miscellaneous.regex_patterns import CommonRegexPatterns

class GlobalChem(object):

    __version__ = "0.9.8"
    __allow_update__ = True

    """

    GlobalChem will be the master class of all variables, as the content store grows we can use this as the parent class.

    """

    def __init__(self):

        pass

    def get_common_regex_patterns(self):

        common_regex_patterns = CommonRegexPatterns()
        return common_regex_patterns.get_common_regex_patterns()

    def get_amino_acids(self):

        amino_acids = AminoAcids()

        amino_acid_smiles = amino_acids.get_amino_acids_smiles()
        amino_acid_smarts = amino_acids.get_amino_acid_smarts()

        return amino_acid_smiles, amino_acid_smarts

    def get_common_organic_solvents(self):

        solvents = CommonOrganicSolvents()

        common_organic_solvents_smiles = solvents.get_organic_solvents_smiles()
        common_organic_solvents_smarts = solvents.get_organic_solvents_smarts()

        return common_organic_solvents_smiles, common_organic_solvents_smarts

    def get_essential_vitamins(self):


        vitamins = Vitamins()

        vitamins_smiles = vitamins.get_vitamin_smiles()
        vitamins_smarts = vitamins.get_vitamin_smarts()

        return vitamins_smiles, vitamins_smarts

    def get_open_smiles_functional_groups(self):


        open_smiles = OpenSmiles()

        functional_groups_smiles = open_smiles.get_open_smiles()
        functional_groups_smarts = open_smiles.get_open_smarts()

        return functional_groups_smiles, functional_groups_smarts

    def get_commonly_used_r_group_replacements(self):

        common_r_groups = CommonRGroupReplacements()

        r_group_smiles = common_r_groups.get_r_group_replacement_smiles()
        r_group_smarts = common_r_groups.get_r_group_replacement_smarts()

        return r_group_smiles, r_group_smarts

    def get_iupac_blue_book_common_functional_groups(self):


        substituents = IUPACBlueBook()
        rings = IUPACBlueBookRings()

        radical_smiles = substituents.get_radical_smiles()
        radical_smarts = substituents.get_radical_smarts()

        rings_smiles = rings.get_rings_smiles()
        rings_smarts = rings.get_rings_smarts()
        
        return radical_smiles, radical_smarts, rings_smiles, rings_smarts

    def get_rings_in_drugs(self):


        rings = RingsInDrugs()

        rings_in_drugs_smiles = rings.get_rings_in_drugs_smiles()
        rings_in_drugs_smarts = rings.get_rings_in_drugs_smarts()

        return rings_in_drugs_smiles, rings_in_drugs_smarts

    def get_common_heterocyclic_rings_phase_2(self):


        rings = Phase2HetereoCyclicRings()

        rings_smiles = rings.get_rings_smiles()
        rings_smarts = rings.get_rings_smarts()

        return rings_smiles, rings_smarts

    def get_common_privileged_scaffolds(self):


        scaffolds = PrivilegedScaffolds()

        privileged_functional_groups_smiles = scaffolds.get_scaffolds_smiles()
        privileged_functional_groups_smarts = scaffolds.get_privileged_smarts()

        return privileged_functional_groups_smiles, privileged_functional_groups_smarts

    def get_common_warhead_covalent_inhibitors(self):

        warheads = CommonWarheadsCovalentInhibitors()

        warhead_smiles = warheads.get_warhead_smiles()
        warhead_smarts = warheads.get_warhead_smarts()

        return warhead_smiles, warhead_smarts

    def get_common_polymer_repeating_units(self):

        monomers = CommonMonomerRepeatingUnits()

        monomers_smiles = monomers.get_monomer_smiles()
        monomers_smarts = monomers.get_monomer_smarts()

        return monomers_smiles, monomers_smarts

    def get_common_electrophilic_warheads_for_kinases(self):

        warheads = ElectrophilicWarheadsForKinases()

        electrophilic_warheads_smiles = warheads.get_warheads_smiles()
        electrophilic_warheads_smarts = warheads.get_warhead_smarts()

        return electrophilic_warheads_smiles, electrophilic_warheads_smarts

    def get_privileged_scaffolds_for_kinase_inhibitors(self):


        scaffolds = PrivilegedKinaseInhibitorScaffolds()

        kinase_inhibitor_smiles = scaffolds.get_scaffolds_smiles()
        kinase_inhibitor_smarts = scaffolds.get_privileged_smarts()

        return kinase_inhibitor_smiles, kinase_inhibitor_smarts

    def get_braf_kinase_inhibitors_for_cancer(self):

        inhibitors = BRAFInhibitorsByPocket()

        ribose_smiles, ribose_smarts = inhibitors.get_ribose_pocket()
        adenine_smiles, adenine_smarts = inhibitors.get_adenine_pocket()

        hydrophobic_smiles, hydrophobic_smarts = inhibitors.get_hydrophobic_pocket()

        type_1_smiles, type_1_smarts = inhibitors.get_type_1_pocket()
        type_2_smiles, type_2_smarts = inhibitors.get_type_2_pocket()

        exposed_to_solvent_smiles, exposed_to_solvent_smarts = inhibitors.get_exposed_to_solvent()

        return (
            ribose_smiles, ribose_smarts,
            adenine_smiles, adenine_smarts,
            hydrophobic_smiles, hydrophobic_smarts,
            type_1_smiles, type_1_smarts,
            type_2_smiles, type_2_smarts,
            exposed_to_solvent_smiles, exposed_to_solvent_smarts
        )

    def get_common_amino_acid_protecting_groups(self):

        protecting_groups = AminoAcidProtectingGroups()

        (
            alpha_amino_removed_by_acid_smiles,
            alpha_amino_removed_by_acid_smarts,
            alpha_amino_removed_by_acid_acronym_smiles,
            alpha_amino_removed_by_acid_acronym_smarts,
            alpha_amino_removed_by_base_smiles,
            alpha_amino_removed_by_base_smarts,
            alpha_amino_removed_by_base_acronym_smiles,
            alpha_amino_removed_by_base_acronym_smarts,
            other_alpha_amino_smiles,
            other_alpha_amino_smarts,
            other_alpha_amino_acronym_smiles,
            other_alpha_amino_acronym_smarts,
        ) =  protecting_groups.get_alpha_amino_acids()

        (
            lys_orn_dap_dab_removed_by_acid_smiles,
            lys_orn_dap_dab_removed_by_acid_smarts,
            lys_orn_dap_dab_removed_by_acid_acronym_smiles,
            lys_orn_dap_dab_removed_by_acid_acronym_smarts,
            lys_orn_dap_dab_removed_by_base_smiles,
            lys_orn_dap_dab_removed_by_base_smarts,
            lys_orn_dap_dab_removed_by_base_acronym_smiles,
            lys_orn_dap_dab_removed_by_base_acronym_smarts,
            other_lys_orn_dap_dab_smiles,
            other_lys_orn_dap_dab_smarts,
            other_lys_orn_dap_dab_acronym_smiles,
            other_lys_orn_dap_dab_acronym_smarts,
        ) = protecting_groups.get_lys_orn_dap_dab_protecting_groups()

        (
            alpha_carboxylic_acid_removed_by_acid_smiles,
            alpha_carboxylic_acid_removed_by_acid_smarts,
            alpha_carboxylic_acid_removed_by_acid_acronym_smiles,
            alpha_carboxylic_acid_removed_by_acid_acronym_smarts,
            alpha_carboxylic_acid_removed_by_base_smiles,
            alpha_carboxylic_acid_removed_by_base_smarts,
            alpha_carboxylic_acid_removed_by_base_acronym_smiles,
            alpha_carboxylic_acid_removed_by_base_acronym_smarts,
            other_alpha_carboxylic_acid_protecting_group_smiles,
            other_alpha_carboxylic_acid_protecting_group_smarts,
            other_alpha_carboxylic_acid_protecting_group_acronym_smiles,
            other_alpha_carboxylic_acid_protecting_group_acronym_smarts,
        ) = protecting_groups.get_alpha_carboxylic_acid_protecting_groups()

        (
            asp_glu_removed_by_acid_smiles,
            asp_glu_removed_by_acid_smarts,
            asp_glu_removed_by_acid_acronym_smiles,
            asp_glu_removed_by_acid_acronym_smarts,
            asp_glu_removed_by_base_smiles,
            asp_glu_removed_by_base_smarts,
            asp_glu_removed_by_base_acronym_smiles,
            asp_glu_removed_by_base_acronym_smarts,
            other_asp_glu_smiles,
            other_asp_glu_smarts,
            other_asp_glu_acronym_smiles,
            other_asp_glu_acronym_smarts,
        ) = protecting_groups.get_asp_glu_protecting_groups()

        (
            amide_backbone_protecting_group_removed_by_acid_smiles,
            amide_backbone_protecting_group_removed_by_acid_smarts,
            amide_backbone_protecting_group_removed_by_acid_acronym_smiles,
            amide_backbone_protecting_group_removed_by_acid_acronym_smarts,
            other_amide_backbone_protecting_group_smiles,
            other_amide_backbone_protecting_group_smarts,
        ) = protecting_groups.get_amide_backbone_protecting_group()

        (
            asn_gln_removed_by_acid_smiles,
            asn_gln_removed_by_acid_smarts,
            asn_gln_removed_by_acid_acronym_smiles,
            asn_gln_removed_by_acid_acronym_smarts,
            arg_gln_removed_by_acid_smiles,
            arg_gln_removed_by_acid_smarts,
            arg_gln_removed_by_acid_acronym_smiles,
            arg_gln_removed_by_acid_acronym_smarts,
            arg_gln_removed_by_base_smiles,
            arg_gln_removed_by_base_smarts,
            arg_gln_removed_by_base_acronym_smiles,
            arg_gln_removed_by_base_acronym_smarts,
            other_arg_gln_smiles,
            other_arg_gln_smarts,
            other_arg_gln_acronym_smiles,
            other_arg_gln_acronym_smarts,
        ) = protecting_groups.get_asn_gln_protecting_groups()

        (
            cys_removed_by_acid_smiles,
            cys_removed_by_acid_smarts,
            cys_removed_by_acid_acronym_smiles,
            cys_removed_by_acid_acronym_smarts,
            cys_removed_by_base_smiles,
            cys_removed_by_base_smarts,
            cys_removed_by_base_acronym_smiles,
            cys_removed_by_base_acronym_smarts,
            other_cys_smiles,
            other_cys_smarts,
            other_cys_acronym_smiles,
            other_cys_acronym_smarts,
        ) = protecting_groups.get_cys_protecting_groups()

        (
            his_removed_by_acid_smiles,
            his_removed_by_acid_smarts,
            his_removed_by_acid_acronym_smiles,
            his_removed_by_acid_acronym_smarts,
            his_removed_by_base_smiles,
            his_removed_by_base_smarts,
            his_removed_by_base_acronym_smiles,
            his_removed_by_base_acronym_smarts,
            other_his_smiles,
            other_his_smarts,
            other_his_acronym_smiles,
            other_his_acronym_smarts,
        ) = protecting_groups.get_his_protecting_groups()

        (
            ser_thr_hyp_removed_by_acid_smiles,
            ser_thr_hyp_removed_by_acid_smarts,
            ser_thr_hyp_removed_by_acid_acronym_smiles,
            ser_thr_hyp_acronym_removed_by_acid_smarts,
            other_ser_thr_hyp_smiles,
            other_ser_thr_hyp_smarts,
            other_ser_thr_hyp_acronym_smiles,
            other_ser_thr_hyp_acronym_smarts,
        ) = protecting_groups.get_ser_thr_hyp_protecting_groups()

        (
            tyr_removed_by_acid_smiles,
            tyr_removed_by_acid_smarts,
            tyr_removed_by_acid_acronym_smiles,
            tyr_removed_by_acid_acronym_smarts,
            other_tyr_protecting_group_smiles,
            other_tyr_protecting_group_smarts,
            other_tyr_protecting_group_acronym_smiles,
            other_tyr_protecting_group_acronym_smarts,
            trp_removed_by_acid_smiles,
            trp_removed_by_acid_smarts,
            trp_removed_by_acid_acronym_smiles,
            trp_removed_by_acid_acronym_smarts,
        ) = protecting_groups.get_tyr_protecting_groups()

        return (
            alpha_amino_removed_by_acid_smiles,
            alpha_amino_removed_by_acid_smarts,
            alpha_amino_removed_by_acid_acronym_smiles,
            alpha_amino_removed_by_acid_acronym_smarts,
            alpha_amino_removed_by_base_smiles,
            alpha_amino_removed_by_base_smarts,
            alpha_amino_removed_by_base_acronym_smiles,
            alpha_amino_removed_by_base_acronym_smarts,
            other_alpha_amino_smiles,
            other_alpha_amino_smarts,
            other_alpha_amino_acronym_smiles,
            other_alpha_amino_acronym_smarts,
            lys_orn_dap_dab_removed_by_acid_smiles,
            lys_orn_dap_dab_removed_by_acid_smarts,
            lys_orn_dap_dab_removed_by_acid_acronym_smiles,
            lys_orn_dap_dab_removed_by_acid_acronym_smarts,
            lys_orn_dap_dab_removed_by_base_smiles,
            lys_orn_dap_dab_removed_by_base_smarts,
            lys_orn_dap_dab_removed_by_base_acronym_smiles,
            lys_orn_dap_dab_removed_by_base_acronym_smarts,
            other_lys_orn_dap_dab_smiles,
            other_lys_orn_dap_dab_smarts,
            other_lys_orn_dap_dab_acronym_smiles,
            other_lys_orn_dap_dab_acronym_smarts,
            alpha_carboxylic_acid_removed_by_acid_smiles,
            alpha_carboxylic_acid_removed_by_acid_smarts,
            alpha_carboxylic_acid_removed_by_acid_acronym_smiles,
            alpha_carboxylic_acid_removed_by_acid_acronym_smarts,
            alpha_carboxylic_acid_removed_by_base_smiles,
            alpha_carboxylic_acid_removed_by_base_smarts,
            alpha_carboxylic_acid_removed_by_base_acronym_smiles,
            alpha_carboxylic_acid_removed_by_base_acronym_smarts,
            other_alpha_carboxylic_acid_protecting_group_smiles,
            other_alpha_carboxylic_acid_protecting_group_smarts,
            other_alpha_carboxylic_acid_protecting_group_acronym_smiles,
            other_alpha_carboxylic_acid_protecting_group_acronym_smarts,
            asp_glu_removed_by_acid_smiles,
            asp_glu_removed_by_acid_smarts,
            asp_glu_removed_by_acid_acronym_smiles,
            asp_glu_removed_by_acid_acronym_smarts,
            asp_glu_removed_by_base_smiles,
            asp_glu_removed_by_base_smarts,
            asp_glu_removed_by_base_acronym_smiles,
            asp_glu_removed_by_base_acronym_smarts,
            other_asp_glu_smiles,
            other_asp_glu_smarts,
            other_asp_glu_acronym_smiles,
            other_asp_glu_acronym_smarts,
            amide_backbone_protecting_group_removed_by_acid_smiles,
            amide_backbone_protecting_group_removed_by_acid_smarts,
            amide_backbone_protecting_group_removed_by_acid_acronym_smiles,
            amide_backbone_protecting_group_removed_by_acid_acronym_smarts,
            other_amide_backbone_protecting_group_smiles,
            other_amide_backbone_protecting_group_smarts,
            asn_gln_removed_by_acid_smiles,
            asn_gln_removed_by_acid_smarts,
            asn_gln_removed_by_acid_acronym_smiles,
            asn_gln_removed_by_acid_acronym_smarts,
            arg_gln_removed_by_acid_smiles,
            arg_gln_removed_by_acid_smarts,
            arg_gln_removed_by_acid_acronym_smiles,
            arg_gln_removed_by_acid_acronym_smarts,
            arg_gln_removed_by_base_smiles,
            arg_gln_removed_by_base_smarts,
            arg_gln_removed_by_base_acronym_smiles,
            arg_gln_removed_by_base_acronym_smarts,
            other_arg_gln_smiles,
            other_arg_gln_smarts,
            other_arg_gln_acronym_smiles,
            other_arg_gln_acronym_smarts,
            cys_removed_by_acid_smiles,
            cys_removed_by_acid_smarts,
            cys_removed_by_acid_acronym_smiles,
            cys_removed_by_acid_acronym_smarts,
            cys_removed_by_base_smiles,
            cys_removed_by_base_smarts,
            cys_removed_by_base_acronym_smiles,
            cys_removed_by_base_acronym_smarts,
            other_cys_smiles,
            other_cys_smarts,
            other_cys_acronym_smiles,
            other_cys_acronym_smarts,
            his_removed_by_acid_smiles,
            his_removed_by_acid_smarts,
            his_removed_by_acid_acronym_smiles,
            his_removed_by_acid_acronym_smarts,
            his_removed_by_base_smiles,
            his_removed_by_base_smarts,
            his_removed_by_base_acronym_smiles,
            his_removed_by_base_acronym_smarts,
            other_his_smiles,
            other_his_smarts,
            other_his_acronym_smiles,
            other_his_acronym_smarts,
            ser_thr_hyp_removed_by_acid_smiles,
            ser_thr_hyp_removed_by_acid_smarts,
            ser_thr_hyp_removed_by_acid_acronym_smiles,
            ser_thr_hyp_acronym_removed_by_acid_smarts,
            other_ser_thr_hyp_smiles,
            other_ser_thr_hyp_smarts,
            other_ser_thr_hyp_acronym_smiles,
            other_ser_thr_hyp_acronym_smarts,
            tyr_removed_by_acid_smiles,
            tyr_removed_by_acid_smarts,
            tyr_removed_by_acid_acronym_smiles,
            tyr_removed_by_acid_acronym_smarts,
            other_tyr_protecting_group_smiles,
            other_tyr_protecting_group_smarts,
            other_tyr_protecting_group_acronym_smiles,
            other_tyr_protecting_group_acronym_smarts,
            trp_removed_by_acid_smiles,
            trp_removed_by_acid_smarts,
            trp_removed_by_acid_acronym_smiles,
            trp_removed_by_acid_acronym_smarts,
        )

    def get_polyfluoroalkyl_substances(self):

        polyfluoroalkyl = EmergingPerFluoroAlkyls()

        polyfluoroalkyl_smiles = polyfluoroalkyl.get_per_smiles()
        polyfluoroalkyl_smarts = polyfluoroalkyl.get_per_smarts()

        return polyfluoroalkyl_smiles, polyfluoroalkyl_smarts

    def get_chemical_adsorption_on_montmorillonite_clays(self):

        adsorption_groups = MontmorilloniteAdsorption()

        functional_groups_smiles = adsorption_groups.get_chemical_smiles()
        functional_groups_smarts = adsorption_groups.get_chemical_smarts()

        return functional_groups_smiles, functional_groups_smarts

       #------------------------- Property Declaration for GlobalChem ---------------------------#
