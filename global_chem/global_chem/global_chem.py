#!/usr/bin/env python3
#
# GlobalChem - Master Object
#
# -----------------------------------

# Base Imports
# ------------

import re
import os
import sys
import os.path
import pprint
from functools import lru_cache

# Reconfigurations
# ----------------

if os.name == 'nt':
    sys.stdin.reconfigure(encoding='utf-8')
    sys.stdout.reconfigure(encoding='utf-8')

# Environment

from global_chem.environment.chemicals_from_biomass import ChemicalsFromBioMass
from global_chem.environment.emerging_perfluoroalkyls import EmergingPerFluoroAlkyls

# Materials

from global_chem.materials.clay.montmorillonite_adsorption import MontmorilloniteAdsorption
from global_chem.materials.polymers.common_monomer_repeating_units import CommonMonomerRepeatingUnits

# Warfare

from global_chem.warfare.organophosphorous_nerve_agents import OrganoPhosphorousNerveAgents

# Education

from global_chem.education.cengage.organic_and_inorganic_bronsted_acids import OrganicAndInorganicBronstedAcids

# Medicinal Chemistry - Cannabis

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

from global_chem.medicinal_chemistry.scaffolds.iupac_blue_book import IUPACBlueBook
from global_chem.medicinal_chemistry.scaffolds.privileged_scaffolds import PrivilegedScaffolds
from global_chem.medicinal_chemistry.scaffolds.common_r_group_replacements import CommonRGroupReplacements

# Proteins Kinases

from global_chem.proteins.kinases.braf.braf_inhibitors import BRAFInhibitors
from global_chem.proteins.kinases.scaffolds.privileged_kinase_inhibitors import PrivilegedKinaseInhibitors

# Organic Synthesis

from global_chem.organic_synthesis.solvents.common_organic_solvents import CommonOrganicSolvents
from global_chem.organic_synthesis.protecting_groups.amino_acid_protecting_groups import AminoAcidProtectingGroups
from global_chem.organic_synthesis.bidendate_phosphine_ligands.nickel_ligands import NickelBidendatePhosphineLigands
# from global_chem.organic_synthesis.named_reactions_in_organic_synthesis.named_reactions_in_organic_synthesis import NamedReactionsInOrganicSynthesis

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

from global_chem.narcotics.pihkal import Pihkal
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

# Miscellaneous

from global_chem.miscellaneous.vitamins import Vitamins
from global_chem.miscellaneous.open_smiles import OpenSmiles
from global_chem.miscellaneous.amino_acids import AminoAcids
from global_chem.miscellaneous.regex_patterns import CommonRegexPatterns

# Sex

from global_chem.sex.exsens.lube import Lube
from global_chem.sex.exsens.exsens_products import ExsensProducts
from global_chem.sex.tainted_sexual_enhancements.tainted_sexual_enhancements import TaintedSexualEnhancements

class Node:

    '''
    Node Object
    '''

    def __init__(self, name, smiles=[], smarts=[], bit_vectors=[], value=None, colour=None):

        self.name = name.split('.')[0]
        self.children = []
        self.parents = []
        self.smiles = smiles
        self.smarts = smarts
        self.bit_vector = bit_vectors
        self.state = [ self.smiles, self.smarts ]
        self.value = None
        self.colour = None

    def add_parent(self, name, smiles=[], smarts=[], bit_vectors=[]):

        '''

        Add the Parent Node

        '''

        self.parents.append(
            Node(name, smiles, smarts, bit_vectors)
        )

    def add_child(self, name, smiles=[], smarts=[], bit_vectors=[]):

        '''

        Add the Child Node

        '''

        self.children.append(
            Node(name, smiles, smarts, bit_vectors)
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

    def set_node_value(self, value):

        '''

        Set the Node Value

        Arguments:
            value: value of the node you would like to set it to | Type: Anything

        '''

        self.value = value

    def get_node_value(self):

        '''

        Get the Node Value

        '''

        return self.value

    # This hack is for the root node to have dummy dictionaries.

    @staticmethod
    def get_smiles():

        smiles = {}

    @staticmethod
    def get_smarts():

        smarts = {}

    @staticmethod
    def get_bit_vector():

        bit_vector = {}

class PrintNode:

    __version__ = '0.0.1'

    '''
    
    Hack for now, 
    
    Print Node function to get the network printed out. Feature first, refactor later. 
    
    '''

    def __init__(self, val, parent=None):

        self.val = val
        self.parent = parent
        self.children = []

        if self.parent:
            self.parent.children.append(self)

    def __str__(self):
        return f"{self.val}"

    def __repr__(self):
        return str(self)

class PrintTreeUtilities(object):

    __version__ = '0.0.1'

    def __init__(self):

        pass

    @staticmethod
    def add_connectors(lst):

        if len(lst) == 1:
            return lst

        startingFanIndex = None
        endingFanIndex = None

        for i in range(0, len(lst)):

            # if we haven't seen an non empty element yet, only add space
            if lst[i].startswith(" ") and startingFanIndex is None:
                continue

            # mark that we just found the first line that doesn't start with indentation
            if startingFanIndex is None:
                startingFanIndex = i

            # update this index as the index of the latest label
            if not lst[i].startswith(" "):
                endingFanIndex = i

        # prepend connectors
        if startingFanIndex == endingFanIndex:
            return lst
        for i in range(len(lst)):
            if i == startingFanIndex:
                lst[i] = f"┌{lst[i]}"
            elif i == endingFanIndex:
                lst[i] = f"└{lst[i]}"
            elif i > startingFanIndex and i < endingFanIndex:
                if lst[i].startswith(" "):
                    lst[i] = f"│{lst[i]}"
                else:
                    lst[i] = f"├{lst[i]}"
            else:
                lst[i] = f" {lst[i]}"
        return lst

    @staticmethod
    def labelFan(node, lst):

        name = str(node)
        if len(lst) == 0:
            lst = [name]
        for i in range(len(lst)):
            # add label to the middle of the fan
            if i == len(lst) // 2:
                lst[i] = f"{name}─{lst[i]}"
            else:
                # push the other labels out with indentation
                indent = " " * len(name)
                lst[i] = f"{indent} {lst[i]}"
        return lst

    @staticmethod
    def connectFans(args):
        union = []
        for n in args:
            union += n
        PrintTreeUtilities.add_connectors(union)
        return union

    @staticmethod
    def printTrees(node):
        def get_repr(node):
            if len(node.children) == 0:
                return [str(node)]
            fans = [get_repr(x) for x in node.children]
            labelled = PrintTreeUtilities.labelFan(node, PrintTreeUtilities.connectFans(fans))
            return labelled
        print(str("\n".join(get_repr(node))))

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

    # NODE CONTRIBUTORS
    # -----------------

    # Thank you to the contributors to the node network that made this project alive.

    __NODES__ = {
        'global_chem': Node,
        'emerging_perfluoroalkyls': EmergingPerFluoroAlkyls,                      # Asuka Orr & Suliman Sharif
        'montmorillonite_adsorption': MontmorilloniteAdsorption,                  # Asuka Orr & Suliman Sharif
        'common_monomer_repeating_units': CommonMonomerRepeatingUnits,            # Suliman Sharif
        'electrophilic_warheads_for_kinases': ElectrophilicWarheadsForKinases,    # Ruibin Liu & Suliman Sharif
        'common_warheads_covalent_inhibitors': CommonWarheadsCovalentInhibitors,  # Shaoqi Zhan & Suliman Sharif
        'rings_in_drugs': RingsInDrugs,                                           # Alexander Mackerell Jr. & Suliman Sharif
        'iupac_blue_book_rings': IUPACBlueBookRings,                              # Suliman Sharif
        'phase_2_hetereocyclic_rings': Phase2HetereoCyclicRings,                  # Suliman Sharif
        'privileged_scaffolds': PrivilegedScaffolds,                              # Suliman Sharif
        'iupac_blue_book': IUPACBlueBook,                                         # Suliman Sharif
        'common_r_group_replacements': CommonRGroupReplacements,                  # Sunhwan Jo & Suliman Sharif
        'braf_inhibitors': BRAFInhibitors,                                        # Aarion Romany & Suliman Sharif
        'privileged_kinase_inhibitors': PrivilegedKinaseInhibitors,               # Suliman Sharif
        'common_organic_solvents': CommonOrganicSolvents,                         # Suliman Sharif
        'amino_acid_protecting_groups': AminoAcidProtectingGroups,                # Aziza Frank & Suliman Sharif
        'schedule_one': ScheduleOne,                                              # Suliman Sharif
        'schedule_two': ScheduleTwo,                                              # Suliman Sharif
        'schedule_three': ScheduleThree,                                          # Suliman Sharif
        'schedule_four': ScheduleFour,                                            # Suliman Sharif
        'schedule_five': ScheduleFive,                                            # Suliman Sharif
        'interstellar_space': InterstellarSpace,                                  # Suliman Sharif
        'vitamins': Vitamins,                                                     # Suliman Sharif
        'open_smiles': OpenSmiles,                                                # Suliman Sharif
        'amino_acids': AminoAcids,                                                # Suliman Sharif
        'pihkal': Pihkal,                                                         # Suliman Sharif
        'nickel_ligands': NickelBidendatePhosphineLigands,                        # Suliman Sharif
        'cimetidine_and_acyclovir': CimetidineAndAcyclovir,                       # Suliman Sharif
        'how_to_live_longer': HowToLiveLonger,                                    # Suliman Sharif
        'monoclonal_antibodies': MonoclonalAntibodies,                            # Asuka Orr & Suliman Sharif
        'lube': Lube,                                                             # Daniel Khavrutskii & Suliman Sharif
        'tainted_sexual_enhancements': TaintedSexualEnhancements,                 # Suliman Sharif
        'salt': Salt,                                                             # Suliman Sharif
        'exsens_products': ExsensProducts,                                        # Rebecca Pinette-Dorin & Suliman Sharif
        'fda_list_one': FDAListOne,                                               # Mike Wostner & Suliman Sharif
        'fda_list_two': FDAListTwo,                                               # Mike Wostner & Suliman Sharif
        'fda_list_three': FDAListThree,                                           # Mike Wostner & Suliman Sharif
        'fda_list_four': FDAListFour,                                             # Mike Wostner & Suliman Sharif
        'fda_list_five': FDAListFive,                                             # Mike Wostner & Suliman Sharif
        'fda_list_six': FDAListSix,                                               # Mike Wostner & Suliman Sharif
        'fda_list_seven': FDAListSeven,                                           # Mike Wostner & Suliman Sharif
        'constituents_of_cannabis_sativa': ConstituentsOfCannabisSativa,          # Ian Jones & Bettina Lier & Suliman Sharif
        'phytocannabinoids': PhytoCannabinoids,                                   # Ian Jones & Bettina Lier & Suliman Sharif
        'organophosphorous_nerve_agents': OrganoPhosphorousNerveAgents,           # Suliman Sharif
        'organic_and_inorganic_bronsted_acids': OrganicAndInorganicBronstedAcids, # Nathaniel McClean & Suliman Sharif
        'chemicals_from_biomass': ChemicalsFromBioMass,                           # Anthony Maiorana & Suliman Sharif
        'common_regex_patterns': CommonRegexPatterns,                             # Chris Burke & Suliman Sharif
    }

    __INCOMPLETE_NODES = {
        # 'named_reactions_in_organic_synthesis': NamedReactionsInOrganicSynthesis # Aziza Frank & Bettina Lier & Suliman Sharif
    }

    def __init__(self, verbose=False):

        self.network = {}
        self.deep_layer_network = {}
        self.root_node = ''
        self.deep_layer_count = 0
        self.verbose = verbose
        self.splitter = '/'

        if os.name == 'nt':
            self.splitter = '\\'

    def check_available_nodes(self):

        '''

        Checks the existing Nodes within the network.

        '''

        return list(self.__NODES__.keys())

    def get_all_nodes(self):

        '''

        Returns all the nodes in no particular order

        '''

        # Skip the head node

        nodes = [ i() for i in list(self.__NODES__.values())[1:] if i().name != 'common_regex_patterns' ]

        return nodes

    def get_node(self, node_key):

        '''
        Return a node based on the key by the user.

        Arguments:
            node_key (String): node keys within the network

        '''

        return self.__NODES__[node_key]

    def get_depth_of_globalchem(self):

        '''

        Returns the Depth of the GlobalChem Tree

        Returns:
            max_depth (Int): Max depth of the object

        '''

        path_objects = []

        if os.name == 'nt':
            absolute_file_path = self.splitter.join(os.path.dirname(os.path.realpath('__file__')).split(self.splitter)[:-1])
        else:
            absolute_file_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])

        for dirpath, dirnames, filenames in os.walk(absolute_file_path):

            for file in filenames:

                if file.endswith('.py') and \
                        '__' not in file and \
                        'cli.py' not in file and \
                        'global_chem.py' not in file:

                    object_path = os.path.join(''.join(dirpath.rsplit(absolute_file_path)), file).split('/')

                    path_objects.append(object_path)

        max_depth = 1

        for i in path_objects:

            depth = len(i) - 1

            if depth > max_depth:
                max_depth = depth

        return max_depth

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
        self.root_node = root_node

        self.network[root_node] = {
            "node_value": Node(
                root_node,
                self.__NODES__[root_node].get_smiles(),
                self.__NODES__[root_node].get_smarts(),
                self.__NODES__[root_node].get_bit_vector()
            ),
            "children": [],
            "parents": [],
            "name": root_node
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
                self.__NODES__[ child_key ].get_bit_vector()
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
                    self.__NODES__[ child_key ].get_smarts(),
                    self.__NODES__[ child_key ].get_bit_vector()
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
                self.__NODES__[ parent_key ].get_smarts(),
                self.__NODES__[ child_key ].get_bit_vector()
            )

        except:

            _ = new_child_node.add_parent(
                parent_key,
            )
        # Step 5

        self.network[ child_key ][ 'parents' ].append(parent_key)

        if self.verbose:
            print (f'Network State Step 5: {self.network}')


    def remove_node(self, name):

        '''
        Rules:
            1.) Cannot Remove a Node that has children.
        Arguments:
            name (String): Name of the node
        '''

        # Check Node Extant

        node = self.network.get(name, None)

        if not node:
            raise GraphNetworkError(
                message=f'Node {name} not found in network',
                errors=''
            )

        # Check the Children

        node_children = node['children']

        if len(node_children) != 0:
            raise GraphNetworkError(
                message=f'Cannot Remove Node {name} because it has children',
                errors=''
            )

        # Remove the Node

        del self.network[name]

        # Remove that navigated to the node

        for key, value in self.network.items():

            value_children = value['children']

            if name in value_children:

                removed_children_items = [ i for i in value_children if i not in [name] ]
                self.network[key]['children'] = removed_children_items

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

    def get_node_bits(self, node_key):

        '''

        Get the Node Bit Vector

        '''

        node = self.network.get(node_key, None)

        if not node:
            raise GraphNetworkError(
                message=f'No Node named {node_key} exists',
                errors=''
            )

        return self.__NODES__[node_key].get_bit_vector()

    def initiate_deep_layer_network(self, root_node='global_chem'):

        '''
        Initiates the Deep Layer Network with the first layer 1 being the root node.
        Arguments:
            root_node (String): Root node to the network.
        '''

        self.root_node = root_node

        self.deep_layer_network[root_node] = {
            "node_value": Node(
                root_node,
                self.__NODES__[root_node].get_smiles(),
                self.__NODES__[root_node].get_smarts()
            ),
            "children": [],
            "parents": [],
            "name": root_node,
            'layer': '1'
        }

        # Increase the count

        self.deep_layer_count += 1

    def add_deep_layer(self, nodes):

        '''

        Arguments:
            nodes (List): Add a Layer of Nodes to the previous parents

        Rules:
            1.) When adding a deep layer add children to all the previous parents.

        Algorithm:
            1.) Fetch all the parents of the current layer
            2.) Add all the node children to the parents.
            3.) Increase the deep layer count.
            4.) All all the parents to the children.

        Graph Algorithm:
            O(cs): Childrens Node
            O(ps): Parents Node
            DN: Deep Network
            DCL: Deep Current Layer

            <F>: Fetch Function
            1.) O(ps) <F> DCL <F> DN
            2.) O(cs) <---- O(ps)
            3.) DCL += 1
            4.) O(ps) ----> O(cs)

        '''

        parents = []

        # Step 1

        for node_key, node_value in self.deep_layer_network.items():

            if self.deep_layer_count == int(node_value['layer']):

                # Step 2

                self.deep_layer_network[node_key]['children'] = nodes
                parents.append(node_key)

        # Step 3

        self.deep_layer_count += 1

        # Step 4

        for child_node in nodes:

            self.deep_layer_network[child_node] = {
                "node_value": Node(
                    child_node,
                    self.__NODES__[child_node].get_smiles(),
                    self.__NODES__[child_node].get_smarts()
                ),
                "children": [],
                "parents": parents,
                "name": child_node,
                'layer': self.deep_layer_count
            }

    def build_global_chem_network(self, print_output=False, debugger=False):

        '''
        Get the GlobalChem Internal Network
        Arguments:
            print_output: Print the output of the globalchem network | Type: Boolean
        '''

        # Initiate the Head

        self.initiate_network()

        # Fetch All the File Paths

        path_objects = []

        if os.name == 'nt':
            absolute_file_path = self.splitter.join(os.path.dirname(os.path.realpath('__file__')).split(self.splitter)[:-1])
        else:
            absolute_file_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])

        for dirpath, dirnames, filenames in os.walk(absolute_file_path):

            for file in filenames:

                if file.endswith('.py') and \
                        '__' not in file and \
                        'cli.py' not in file and \
                        'global_chem.py' not in file:

                    object_path = os.path.join(''.join(dirpath.rsplit(absolute_file_path)), file).split(self.splitter)

                    if debugger:
                        print ("Node Object Paths: %s: " % object_path)

                    path_objects.append(object_path)
        
        # Add the objects recursively

        for chemical_object in path_objects:

            chemical_object.reverse()
            chemical_object = chemical_object + ['global_chem']

            while len(chemical_object) > 0:

                parent = chemical_object.pop().split('.')[0]
                previous_child = parent

                if len(chemical_object) != 0:
                    child = chemical_object[-1].split('.')[0]
                    self.add_node(parent, child)
                else:
                    self.add_node(previous_child, parent)

        self.network[ 'global_chem' ][ 'children' ] = list(set(self.network[ 'global_chem' ][ 'children' ]))

        # Pretty Print the Objects

        if print_output:
            pretty_printer = pprint.PrettyPrinter()
            pretty_printer.pprint(self.network)

    def get_all_names(self):

        '''
        Fetches all the names in the network
        '''

        names = []

        for node_key, node_value in self.__NODES__.items():

            if node_key != 'global_chem' and node_key != 'common_regex_patterns':

                names.append(list(node_value.get_smiles().keys()))

        names = sum(names, [])

        return names

    def get_all_smiles(self):

        '''
        Fetches all the smiles in the network
        Returns:
            smiles (List): Full list of SMILES in the network.
        '''

        smiles = []

        for node_key, node_value in self.__NODES__.items():

            if node_key != 'global_chem' and node_key != 'common_regex_patterns':

                smiles.append(list(node_value.get_smiles().values()))

        smiles = sum(smiles, [])

        return smiles

    def get_all_smarts(self):

        '''
        Fetches all the smarts in the network
        Returns:
            smarts (List): Full list of SMARTS in the network.
        '''

        smarts = []

        for node_key, node_value in self.__NODES__.items():

            if node_key != 'global_chem' and node_key != 'common_regex_patterns':

                smarts.append(list(node_value.get_smarts().values()))

        smarts = sum(smarts, [])

        return smarts

    def get_all_bits(self):

        '''

        Fetches all the smarts in the network

        Returns:
            bits (List): Full list of SMARTS in the network.

        '''

        bits = []

        for node_key, node_value in self.__NODES__.items():

            if node_key != 'global_chem' and node_key != 'common_regex_patterns':

                bits.append(list(node_value.get_bit_vector().values()))

        bits = sum(bits, [])

        return bits

    def get_smiles_by_iupac(
            self,
            iupac_key,
            distance_tolerance=4,
            return_partial_definitions=False,
            reconstruct_smiles=False,
        ):

        '''

        Get the SMILES by the IUPAC.

        Arguments:
            iupac_key (String): Key for the iupac.
            distance_tolerance (Int): Distance tolerance for Levenshetin Distances.
            return_partial_definitions (Bool): Return the Levenshetin Distances between words
            reconstruct_smiles (Bool): Reconstruct the Chemical Space SMILES from the IUPAC
        Returns:
            definition (Dict): A definition of the IUPAC.

        '''

        if reconstruct_smiles:

            reconstructed_smiles = ''

            iupac_keys = iupac_key.lower().split('-')
            iupac_keys = [re.sub(r'\b[^a-zA-Z]\b', '', iupac_key) for iupac_key in iupac_keys]
            iupac_keys = [iupac_key for iupac_key in iupac_keys if not any(char.isdigit() for char in iupac_key) and iupac_key]


            for node_key, node_value in self.__NODES__.items():

               if node_key == 'global_chem' or node_key == 'common_regex_patterns':
                   continue

               entity = node_value()
               names = entity.get_smiles()

               for name, smiles in names.items():

                   # Sanitize the Node Key

                   name = name.lower()
                   name = re.sub(r'[^a-zA-Z]', '', name)

                   for iupac_key in iupac_keys:

                       distance = self.levenshtein_distance(iupac_key, name)

                       if 0 <= distance <= distance_tolerance:

                           reconstructed_smiles += smiles + '.'

            return reconstructed_smiles[:-1]
        else:

            iupac_key = iupac_key.lower()
            iupac_key = re.sub(r'[^a-zA-Z]', '', iupac_key)

            definitions = []
            exact_definition = {}

            for node_key, node_value in self.__NODES__.items():

                if node_key == 'global_chem' or node_key == 'common_regex_patterns':
                    continue

                network_path = node_value.__module__
                entity = node_value()
                names = entity.get_smiles()

                for name, smiles in names.items():

                    definition = {}

                    # Sanitize the Node Key

                    name = name.lower()
                    name = re.sub(r'[^a-zA-Z]', '', name)

                    distance = self.levenshtein_distance(iupac_key, name)

                    if 0 <= distance <= distance_tolerance:

                        definition[name] = smiles
                        definition['network_path'] = network_path
                        definition['levenshtein_distance'] = distance

                        if distance == 0:
                            exact_definition = definition

                        definitions.append(definition)

            if return_partial_definitions:
                return definitions
            else:
                return exact_definition

    def set_node_value(self, node_key, value):

        '''
        Set the Node value
        '''

        node = self.network.get(node_key, None)

        if not node:
            raise GraphNetworkError(
                message=f'No Node named {node_key} exists',
                errors=''
            )

        self.network[node_key]['node_value'].set_node_value(value)

    def get_node_value(self, node_key):

        '''
        Get the Node Value
        Arguments:
            node_key (value): Get the Node Key.
        '''

        node = self.network.get(node_key, None)

        if not node:
            raise GraphNetworkError(
                message=f'No Node named {node_key} exists',
                errors=''
            )

        return self.network[node_key]['node_value'].get_node_value()

    def compute_common_score(self, iupac_name, verbose=False):

        '''
        Compute the Common score for an IUPAC name
        Arguments:
            iupac_name (String): iupac name in question
            verbose (Boolean): verbose flag for the common score
        Returns:
            common_score (Float): Common Score computed by GlobalChem Algorithm
        Common Score Algorithm:
            1.) Data mine the current state of GlobalChem
            2.) Get the Object Weights of Each mention
            3.) Determine the Mention Weight
            4.) Sum the Weights and That's How common it is.
        '''

        # Step 1

        iupac_names = []
        object_names =[]

        for node_key, node_value in self.__NODES__.items():

            if node_key != 'global_chem' and node_key != 'common_regex_patterns':

                object_names.append(node_key)
                iupac_names.append(list(node_value.get_smiles().keys()))

        # Step 2

        total_names = len(sum(iupac_names, []))
        object_weights = [ len(i) / total_names for i in iupac_names]

        # Step 3

        mentioned_weights = []
        total_weights = []

        for i in range(0, len(iupac_names)):

            names = iupac_names[i]

            mentioned = 0

            for name in names:

                if iupac_name in name:

                    mentioned += 1

            mentioned_weight = mentioned * object_weights[i]
            mentioned_weights.append(mentioned_weight)

            # print ("Object Name: %s" % object_names[i])
            # print ("Mentioned Weight: %s" % mentioned_weight)

        # Step 4

        common_score = sum(mentioned_weights)

        if verbose:
            print ("GlobalChem Common Score: %s" % common_score)

        return common_score

    def print_globalchem_network(self):


        '''
        Outputs the Network into the terminal into a nice format.
        '''

        # Initialize Print node network

        _PRINT_NODE_KEY = {}

        path_objects = []

        if os.name == 'nt':
            absolute_file_path = self.splitter.join(os.path.dirname(os.path.realpath('__file__')).split(self.splitter)[:-1])
        else:
            absolute_file_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])

        for dirpath, dirnames, filenames in os.walk(absolute_file_path):

            for file in filenames:

                if file.endswith('.py') and \
                        '__' not in file and \
                        'cli.py' not in file and \
                        'global_chem.py' not in file:

                    object_path = os.path.join(''.join(dirpath.rsplit(absolute_file_path)), file).split(self.splitter)
                    path_objects.append(object_path)

        # Add the objects recursively

        for chemical_object in path_objects:

            chemical_object.reverse()
            chemical_object = chemical_object + ['global_chem']

            while len(chemical_object) > 0:

                parent = chemical_object.pop().split('.')[0]
                previous_child = parent

                if len(chemical_object) != 0:

                    child = chemical_object[-1].split('.')[0]

                    if parent not in _PRINT_NODE_KEY:

                        _PRINT_NODE_KEY[ parent ] = PrintNode(parent)

                    if child not in _PRINT_NODE_KEY:

                        _PRINT_NODE_KEY[ child ] = PrintNode( child, _PRINT_NODE_KEY[ parent ] )

                else:

                    if not previous_child:

                        _PRINT_NODE_KEY[ previous_child ] = PrintNode(previous_child)

                    if not parent:

                        _PRINT_NODE_KEY[ parent ] = PrintNode(child, _PRINT_NODE_KEY[parent])

        print(str(PrintTreeUtilities.printTrees(_PRINT_NODE_KEY['global_chem'])))

    def print_deep_network(self):


        '''
        Outputs the Deep Network into the terminal into a nice format.
        '''


        _DEEP_NETWORK_KEY = {}

        _DEEP_NETWORK_KEY[ self.root_node ] = PrintNode(self.root_node)

        previous_parents = []

        for i in range(self.deep_layer_count):

            layer_number = i + 2

            for node_key, node_value in self.deep_layer_network.items():

                if node_value['layer'] == layer_number:

                    if layer_number == 2:

                        _DEEP_NETWORK_KEY[ node_key ] = PrintNode(node_key, _DEEP_NETWORK_KEY[ self.root_node ])

                        previous_parents.append(node_key)

                    else:

                        for previous_parent in previous_parents:

                            _DEEP_NETWORK_KEY[ node_key ] = PrintNode(node_key, _DEEP_NETWORK_KEY[ previous_parent ])

        print(PrintTreeUtilities.printTrees(_DEEP_NETWORK_KEY[ self.root_node ]))

    def to_tsv(self, out_file = 'global_chem.tsv'):


        '''
        Convert the network into a TSV file in the format:
            IUPAC, SMILES, NODE NAME, TREE PATH
        '''

        master_category_keys = {
            'organic_chemistry': ['narcotics', 'organic_synthesis', 'medicinal_chemistry', 'proteins', 'miscellaneous'],
            'environmental_chemistry': ['environment', 'interstellar_space'],
            'materials_chemistry': ['materials'],
            'pharmaceutical_sciences': ['formulation']
        }

        out_file = open(out_file, 'w')

        all_nodes = self.get_all_nodes()

        for node in all_nodes:

            node_name = node.name
            tree_path = node.__module__

            category = ''
            sub_category = tree_path.split('.')[1]

            for key, value in master_category_keys.items():
                if sub_category in value:
                    category = key

            for key, value in node.get_smiles().items():

                out_file.write(f'{key}\t{value}\t{node_name}\t{category}\t{tree_path}\n')

        out_file.close()

    @staticmethod
    def levenshtein_distance(iupac_a, iupac_b):

        '''

        This function will calculate the levenshtein distance between two IUPAC strings

        params:
            a (String) : The first string you want to compare
            b (String) : The second string you want to compare

        returns:
            This function will return the distnace between string a and b.

        '''

        @lru_cache(None)
        def min_dist(string_1, string_2):

            if string_1 == len(iupac_a) or string_2 == len(iupac_b):
                return len(iupac_a) - string_1 + len(iupac_b) - string_2

            if iupac_a[string_1] == iupac_b[string_2]:
                return min_dist(string_1 + 1, string_2 + 1)

            return 1 + min(
                min_dist(string_1, string_2 + 1),
                min_dist(string_1 + 1, string_2),
                min_dist(string_1 + 1, string_2 + 1),
            )

        return min_dist(0, 0)