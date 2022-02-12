#!/usr/bin/env python3
#
# Intialization of the package
#
# ----------------------------


# Environment

from global_chem.environment.emerging_perfluoroalkyls import EmergingPerFluoroAlkyls

# Materials

from global_chem.materials.clay.montmorillonite_adsorption import MontmorilloniteAdsorption
from global_chem.materials.polymers.common_monomer_repeating_units import CommonMonomerRepeatingUnits

# Medicinal Chemistry - Warheads

from global_chem.medicinal_chemistry.warheads.electrophillic_warheads_for_kinases import ElectrophilicWarheadsForKinases
from global_chem.medicinal_chemistry.warheads.common_warheads_covalent_inhibitors import CommonWarheadsCovalentInhibitors

# Medicinal Chemistry - Rings

from global_chem.medicinal_chemistry.rings.rings_in_drugs import RingsInDrugs
from global_chem.medicinal_chemistry.rings.iupac_blue_book_rings import IUPACBlueBookRings
from global_chem.medicinal_chemistry.rings.phase_2_hetereocyclic_rings import Phase2HetereoCyclicRings

# Medicinal Chemistry - Scaffolds

from global_chem.medicinal_chemistry.scaffolds.privileged_scaffolds import PrivilegedScaffolds
from global_chem.medicinal_chemistry.scaffolds.iupac_blue_book_substituents import IUPACBlueBook
from global_chem.medicinal_chemistry.scaffolds.common_r_group_replacements import CommonRGroupReplacements

# Proteins Kinases

from global_chem.proteins.kinases.braf.inhibitors import BRAFInhibitors
from global_chem.proteins.kinases.scaffolds.privileged_kinase_inhibtors import PrivilegedKinaseInhibitors

# Organic Synthesis

from global_chem.organic_synthesis.solvents.common_organic_solvents import CommonOrganicSolvents
from global_chem.organic_synthesis.protecting_groups.amino_acid_protecting_groups import AminoAcidProtectingGroups

# Narcotics

from global_chem.narcotics.schedule_one import ScheduleOne
from global_chem.narcotics.schedule_two import ScheduleTwo
from global_chem.narcotics.schedule_three import ScheduleThree
from global_chem.narcotics.schedule_four import ScheduleFour
from global_chem.narcotics.schedule_five import ScheduleFive

# Interstellar Space

from global_chem.interstellar_space.interstellar_space import InterstellarSpace

# Miscellaneous

from global_chem.miscellaneous.vitamins import Vitamins
from global_chem.miscellaneous.open_smiles import OpenSmiles
from global_chem.miscellaneous.amino_acids import AminoAcids
from global_chem.miscellaneous.regex_patterns import CommonRegexPatterns

from global_chem.global_chem import GlobalChem

__all__ = ['GlobalChem']

__name__ = 'GlobalChem'