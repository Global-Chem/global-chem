#!/usr/bin/env python3
#
# GlobalChem - Master Object
#
# -----------------------------------

# Base Imports

import os
import os.path
import pprint

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

class Node:

    '''

    Node Object

    '''

    def __init__(self, name, smiles=[], smarts=[]):

        self.name = name.split('.')[0]
        self.children = []
        self.parents = []
        self.smiles = smiles
        self.smarts = smarts
        self.state = [ self.smiles, self.smarts ]

    def add_parent(self, name, smiles=[], smarts=[]):

        '''

        Add the Parent Node

        '''

        self.parents.append(
            Node(name, smiles, smarts)
        )

    def add_child(self, name, smiles=[], smarts=[]):

        '''

        Add the Child Node

        '''

        self.children.append(
            Node(name, smiles, smarts)
        )

    def get_node_state(self):

        '''

        Get the Node State

        Returns:
            state (Node): Node state

        '''

        return self.state

    def set_node_state(self):

        '''

        Set the Node State with the children and the parents.

        '''

        self.state = [ self.children, self.parents ]

    def print_stat(self):

        """

        Print statistics of the node

        """
        print("Children: %s" % self.children)
        print("Parents: %s" % self.parents)

    # This hack is for the root node to have dummy dictionaries.

    @staticmethod
    def get_smiles():

        smiles = {}

    @staticmethod
    def get_smarts():

        smarts = {}


class GraphNetworkError(Exception):

    __version_error_parser__ = "1.0.4"
    __allow_update__ = False

    '''
    
    Raise the Network Error if the Node cannot be found.
    
    '''
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class GlobalChem(object):

    __version__ = "1.0.4"
    __allow_update__ = False

    """

    GlobalChem will be the master class of all variables, as the content store grows we can use this as the parent class.

    """

    __NODES__ = {
        'global_chem': Node,
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

    def __init__(self, verbose=False):

        self.network = {}
        self.verbose = verbose

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

    def initiate_network(self, root_node='global_chem'):

        '''

        Initiates a Tree Node Network with Root Node

        Arguments:
            root_node: Add the root node to the network | Type: String

        Objects Initialized:
            network: network object

        '''

        self.network = {}

        self.network[root_node] = {
            "node_value": Node(
                root_node,
                self.__NODES__[root_node].get_smiles(),
                self.__NODES__[root_node].get_smarts()
            ),
            "children": [],
            "parents": [],
            "name": root_node.split('.')[0]
        }

    def add_node(self, parent_key, child_key):

        '''

        Add a node into the network

        Rules:
            1.) There must be a connection from one node to the next. Nodes in space cannot exist.

        Algorithm:

            1.) First Find the Node in the Network
            2.) Add the Child to the Parent
            3.) Add the Child to the Network
            4.) Add the Parent to the Child
            5.) Add the Parent to the Network

        Graph Algorithm:

            O(c): Child Node
            O(p): Parent Node
            N: Network

            1.) O(p) <---- N
            2.) O(c) <---- O(p)
            3.) O(c) <---- N
            4.) O(p) ----> O(c)
            5.) O(p) ----> N

        Errors:

            GraphNetworkError Step 1: If the Parent Node was never found
            GraphNetworkError Step 4: If the Newly Child node was never added.

        '''

        # Step 1

        try:
            parent_node = self.network[ parent_key ][ 'node_value' ]
        except:
            raise GraphNetworkError(
                message=f'Step 1 Node: {parent_key} Not Found',
                errors=''
            )

        # Step 2

        try:
            _ = parent_node.add_child(
                child_key,
                self.__NODES__[ child_key ].get_smiles(),
                self.__NODES__[ child_key ].get_smarts(),
            )
        except:
            _ = parent_node.add_child(
                child_key,
            )

        # Step 3

        self.network[ parent_key ][ 'children' ].append(child_key)

        try:
            self.network[ child_key ] = {
                "node_value": Node(
                    child_key,
                    self.__NODES__[ child_key ].get_smiles(),
                    self.__NODES__[ child_key ].get_smarts()
                ),
                "children": [],
                "parents": [],
                "name": child_key.split('.')[0]
            }

        except:

            self.network[ child_key ] = {
                "node_value": Node(
                    child_key,
                    {},
                    {}
                ),
                "children": [],
                "parents": [],
                "name": child_key.split('.')[0]
            }

        if self.verbose:

            print (f'Network State Step 2:')
            print ('----------------------')

            for node_key, node_value in self.network.items():
                print ('Node')
                print (f'{node_key}: {node_value}')
                print ('\n')

        # Step 4

        try:
            new_child_node = self.network[ child_key ]['node_value']

        except:
            raise GraphNetworkError(
                message=f'Node: {child_key} Not Added Correctly',
                errors=''
            )

        try:
            _ = new_child_node.add_parent(
                parent_key,
                self.__NODES__[ parent_key ].get_smiles(),
                self.__NODES__[ parent_key ].get_smarts()
            )

        except:

            _ = new_child_node.add_parent(
                parent_key,
            )
        # Step 5

        self.network[ child_key ][ 'parents' ].append(parent_key)

        if self.verbose:
            print (f'Network State Step 5: {self.network}')


    def get_node(self, name):

        '''

        Get the Node

        Arguments:
            name: Parent key of the node | Type: String

        Returns:
            node: Returns the node.

        Errors:
            GraphNetworkError: If the node doesn't exist.

        '''

        # Hack for the python file: Need to fix later on.

        python_name = name + '.py'

        try:
            node = self.network[name]
        except:
            try:
                node = self.network[python_name]
            except:
                raise GraphNetworkError(
                    message=f'No Node named {name} exists',
                    errors=''
                )

        return node

    def get_node_smiles(self, node_key):

        '''

        Get the Node SMILES

        Arguments:
            node_key: Node key to query | Type: String

        Returns:
            smiles_dict: Dictionary of the SMILES for that object.

        '''

        node = self.network.get(node_key, None)

        if not node:
            raise GraphNetworkError(
                message=f'No Node named {node_key} exists',
                errors=''
            )

        return self.__NODES__[node_key].get_smiles()

    def get_node_smarts(self, node_key):

        '''

        Get the Node SMARTS

        Arguments:
            node_key: Node key to query | Type: String

        Returns:
            smarts_dict: Dictionary of the SMARTS for that object.

        '''

        node = self.network.get(node_key, None)

        if not node:
            raise GraphNetworkError(
                message=f'No Node named {node_key} exists',
                errors=''
            )

        return self.__NODES__[node_key].get_smarts()

    def build_global_chem_network(self, print_output=False):

        '''

        Get the GlobalChem Internal Network

        Arguments:
            print_output: Print the output of the globalchem network | Type: Boolean

        '''

        # Initiate the Head

        self.initiate_network()

        # Fetch All the File Paths

        path_objects = []

        for dirpath, dirnames, filenames in os.walk("."):

            for file in filenames:

                if file.endswith('.py') and \
                        '__' not in file and \
                        'cli.py' not in file and \
                        'global_chem.py' not in file:

                    object_path = os.path.join(dirpath, file).split('/')[1:]
                    path_objects.append(object_path)

        # Add the objects recursively

        for chemical_object in path_objects:

            chemical_object.reverse()
            chemical_object = chemical_object + ['global_chem']

            while len(chemical_object) > 0:

                parent = chemical_object.pop()
                previous_child = parent

                if len(chemical_object) != 0:
                    child = chemical_object[-1]
                    self.add_node(parent, child)
                else:
                    self.add_node(previous_child, parent)

        self.network[ 'global_chem' ][ 'children' ] = list(set(self.network[ 'global_chem' ][ 'children' ]))

        # Pretty Print the Objects

        pretty_printer = pprint.PrettyPrinter()

        pretty_printer.pprint(self.network)
