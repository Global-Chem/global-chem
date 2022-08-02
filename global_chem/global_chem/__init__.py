#!/usr/bin/env python3
#
# Intialization of the package
#
# ----------------------------


# Environment

from global_chem.environment.chemicals_from_biomass import ChemicalsFromBioMass
from global_chem.environment.emerging_perfluoroalkyls import EmergingPerFluoroAlkyls

# Materials

from global_chem.materials.clay.montmorillonite_adsorption import MontmorilloniteAdsorption
from global_chem.materials.polymers.common_monomer_repeating_units import CommonMonomerRepeatingUnits

# Education

from global_chem.education.cengage.organic_and_inorganic_bronsted_acids import OrganicAndInorganicBronstedAcids

# Warfare

from global_chem.warfare.organophosphorous_nerve_agents import OrganoPhosphorousNerveAgents

# Medicinal Chemistry - Cannabinoids

from global_chem.medicinal_chemistry.cannabinoids.phytocannabinoids import PhytoCannabinoids
from global_chem.medicinal_chemistry.cannabinoids.constituents_of_cannabis_sativa import ConstituentsOfCannabisSativa

# Medicinal Chemistry - International

from global_chem.medicinal_chemistry.chinese.how_to_live_longer import HowToLiveLonger

# Medicinal Chemistry - Warheads

from global_chem.medicinal_chemistry.warheads.electrophilic_warheads_for_kinases import ElectrophilicWarheadsForKinases
from global_chem.medicinal_chemistry.warheads.common_warheads_covalent_inhibitors import CommonWarheadsCovalentInhibitors

# Medicinal Chemistry - Rings

from global_chem.medicinal_chemistry.rings.rings_in_drugs import RingsInDrugs
from global_chem.medicinal_chemistry.rings.iupac_blue_book_rings import IUPACBlueBookRings
from global_chem.medicinal_chemistry.rings.phase_2_hetereocyclic_rings import Phase2HetereoCyclicRings

# Medicinal Chemistry - Scaffolds

from global_chem.medicinal_chemistry.scaffolds.privileged_scaffolds import PrivilegedScaffolds
from global_chem.medicinal_chemistry.scaffolds.iupac_blue_book import IUPACBlueBook
from global_chem.medicinal_chemistry.scaffolds.common_r_group_replacements import CommonRGroupReplacements

# Proteins Kinases

from global_chem.proteins.kinases.braf.braf_inhibitors import BRAFInhibitors
from global_chem.proteins.kinases.scaffolds.privileged_kinase_inhibitors import PrivilegedKinaseInhibitors

# Organic Synthesis

from global_chem.organic_synthesis.solvents.common_organic_solvents import CommonOrganicSolvents
from global_chem.organic_synthesis.protecting_groups.amino_acid_protecting_groups import AminoAcidProtectingGroups
from global_chem.organic_synthesis.bidendate_phosphine_ligands.nickel_ligands import NickelBidendatePhosphineLigands

# Food

from global_chem.food.salt.salt import Salt
from global_chem.food.color_additives.fda_list_one import FDAListOne
from global_chem.food.color_additives.fda_list_two import FDAListTwo
from global_chem.food.color_additives.fda_list_three import FDAListThree
from global_chem.food.color_additives.fda_list_four import FDAListFour
from global_chem.food.color_additives.fda_list_five import FDAListFive
from global_chem.food.color_additives.fda_list_six import FDAListSix
from global_chem.food.color_additives.fda_list_seven import FDAListSeven


# Narcotics

from global_chem.narcotics.schedule_one import ScheduleOne
from global_chem.narcotics.schedule_two import ScheduleTwo
from global_chem.narcotics.schedule_three import ScheduleThree
from global_chem.narcotics.schedule_four import ScheduleFour
from global_chem.narcotics.schedule_five import ScheduleFive

# Interstellar Space

from global_chem.interstellar_space.interstellar_space import InterstellarSpace

# Biopharmaceutics - Excipients

from global_chem.formulation.excipients.monoclonal_antibodies.monoclonal_antibodies import MonoclonalAntibodies
from global_chem.formulation.excipients.biopharmaceutics_class_three.cimetidine_and_acyclovir import CimetidineAndAcyclovir

# Sex

from global_chem.sex.exsens.lube import Lube
from global_chem.sex.exsens.exsens_products import ExsensProducts
from global_chem.sex.tainted_sexual_enhancements.tainted_sexual_enhancements import TaintedSexualEnhancements

# Miscellaneous

from global_chem.miscellaneous.vitamins import Vitamins
from global_chem.miscellaneous.open_smiles import OpenSmiles
from global_chem.miscellaneous.amino_acids import AminoAcids
from global_chem.miscellaneous.regex_patterns import CommonRegexPatterns

from global_chem.global_chem import GlobalChem

__all__ = ['GlobalChem']

__name__ = 'GlobalChem'
