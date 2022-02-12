#!/usr/bin/env python3
#
# GlobalChem - Content Variable Store
#
# -----------------------------------

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

# InterstellarSpace

from global_chem.interstellar_space.interstellar_space import InterstellarSpace

# Miscellaneous

from global_chem.miscellaneous.vitamins import Vitamins
from global_chem.miscellaneous.open_smiles import OpenSmiles
from global_chem.miscellaneous.amino_acids import AminoAcids
from global_chem.miscellaneous.regex_patterns import CommonRegexPatterns

class Node:

    '''

    Node Object

    '''

    def __init__(self, value, internal_object):

        self.children = []
        self.value = value
        self.internal_object = internal_object

    def add_child(self, value):
        self.children.append(Node(value, self.internal_object))

    def __repr__(self):
        classname = type(self).__name__

        return (
            f'{classname}({self.value!r}, {self.children})' if self.children else  f'{classname}({self.value!r})')

    def print_stat(self):

        """

        Print statistics of the node

        """
        print("Children: %s" % self.children)
        print("Values: %s" % self.value)

class GlobalChem(object):

    __version__ = "1.0.4"
    __allow_update__ = False

    """

    GlobalChem will be the master class of all variables, as the content store grows we can use this as the parent class.

    """

    __NODES__ = {
        'emerging_perfluoro_alkyls': EmergingPerFluoroAlkyls,
        'montmorillonite_adsorption': MontmorilloniteAdsorption,
        'common_monomer_repeating_units': CommonMonomerRepeatingUnits,
        'electrophilic_warheads_for_kinases': ElectrophilicWarheadsForKinases,
        'common_warhead_covalent_inhibitors': CommonWarheadsCovalentInhibitors,
        'rings_in_drugs': RingsInDrugs,
        'iupac_blue_book_rings': IUPACBlueBookRings,
        'phase_2_hetereocyclic_rings': Phase2HetereoCyclicRings,
        'privileged_scaffolds': PrivilegedScaffolds,
        'iupac_blue_book': IUPACBlueBook,
        'common_rgroup_replacements': CommonRGroupReplacements,
        'braf_inhibitors': BRAFInhibitors,
        'privileged_kinase_inhibitor_scaffolds': PrivilegedKinaseInhibitors,
        'common_organic_solvents': CommonOrganicSolvents,
        'amino_acid_protecting_groups': AminoAcidProtectingGroups,
        'schedule_one': ScheduleOne,
        'schedule_two': ScheduleTwo,
        'schedule_three': ScheduleThree,
        'schedule_four': ScheduleFour,
        'schedule_five': ScheduleFive,
        'interstellar_space': InterstellarSpace,
        'vitamins': Vitamins,
        'open_smiles': OpenSmiles,
        'amino_acids': AminoAcids,
        'common_regex_patterns': CommonRegexPatterns,
    }

    def __init__(self):

        self.network = {}

        pass

    def check_available_nodes(self):

        '''

        Checks the existing Nodes within the network.

        '''

        return list(self.__NODES__.keys())

    def get_all_nodes(self):

        '''

        Returns all the nodes in no particular order

        '''

        return list(self.__NODES__.values())

    def get_node(self, node_key):

        '''

        Return a node based on the key by the user.

        Arguments:

            node_key (String): node keys within the network
        '''

        return self.__NODES__[node_key]

    def get_nodes(self, node_keys):

        '''

        Arguments:
            node_keys (List): return the list of nodes based on the user keys.

        Returns
            nodes (List): List of the node objects

        '''

        nodes = []

        for key in node_keys:

            nodes.append(self.__NODES__[key])

        return nodes

    def initiate_network(self):

        '''

        Initiates a Tree Node Network

        Objects Initialized:
            network: network object

        '''

        self.network = {}

        # Place holder for now for new features.

    def add_node(self, node_key, parent=None):

        '''

        Add a node into the network

        '''

        if not parent:
            self.network[node_key] = Node(node_key, self.__NODES__[node_key])

        else:
            branch = self.network.get(parent, None)
            if not branch:
                raise KeyError(f'No Node named {parent} exists')

            branch.add_child(node_key)

    def get_parent(self, parent):

        '''

        Get the Root node of the parent.

        '''

        node_key = self.network.get(parent, None)

        if not node_key:
            raise KeyError(f'No Node named {parent} exists')

        return node_key.children
