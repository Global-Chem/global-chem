#!/usr/bin/env python3
#
# Intialization of the package
#
# ----------------------------

from global_chem.global_chem import GlobalChem

# Environment

from global_chem.environment.emerging_perfluoroalkyls import EmergingPerFluoroAlkyls

# Medicinal Chemistry

from global_chem.medicinal_chemistry.scaffolds.privileged_scaffolds import PrivilegedScaffolds
from global_chem.medicinal_chemistry.scaffolds.iupac_blue_book_substituents import IUPACBlueBook
from global_chem.medicinal_chemistry.scaffolds.common_r_group_replacements import CommonRGroupReplacements

from global_chem.medicinal_chemistry.warheads.electrophillic_warheads_for_kinases import ElectrophilicWarheadsForKinases
from global_chem.medicinal_chemistry.warheads.common_warheads_covalent_inhibitors import CommonWarheadsCovalentInhibitors

from global_chem.medicinal_chemistry.rings.rings_in_drugs import RingsInDrugs
from global_chem.medicinal_chemistry.rings.phase_2_hetereocyclic_rings import Phase2HetereoCyclicRings
from global_chem.medicinal_chemistry.rings.iupac_blue_book_rings import IUPACBlueBookRings

# Organic Synthesis

from global_chem.organic_synthesis.protecting_groups.amino_acid_protecting_groups import AminoAcidProtectingGroups
from global_chem.organic_synthesis.solvents.common_organic_solvents import CommonOrganicSolvents

# Proteins

from global_chem.proteins.kinases.scaffolds.privileged_scaffolds import PrivilegedKinaseInhibitorScaffolds
from global_chem.proteins.kinases.braf.common_groups_by_pocket import BRAFInhibitorsByPocket

# Materials

from global_chem.materials.polymers.common_monomer_repeating_units import CommonMonomerRepeatingUnits
from global_chem.materials.clay.montmorillonite_adsorption import MontmorilloniteAdsorption

# Miscellaneous

from global_chem.miscellaneous.vitamins import Vitamins
from global_chem.miscellaneous.amino_acids import AminoAcids
from global_chem.miscellaneous.open_smiles import OpenSmiles
from global_chem.miscellaneous.regex_patterns import CommonRegexPatterns

name = 'GlobalChem'