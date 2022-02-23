#!/usr/bin/env python3
#
# GlobalChemExtensions - Networkx Adapter
#
# ---------------------------------------

# Imports
# -------

import networkx as nx

# Graphing Imports
# ----------------

import matplotlib.pyplot as plt


class NetworkxAdapter(object):

    __version__ = '0.0.1'

    def __init__(self):

        self.network = {}
        self.networkx_graph = nx.Graph()

    def convert(self, network):

        '''

        Convert to a networkx object

        '''

        self.network = network

        # First Establish all the nodes

        for node_key, node_value in self.network.items():

            self.networkx_graph.add_node(node_value['name'])

        # Second Add Child Connections

        for node_key, node_value in self.network.items():

            parent_name = node_value['name']
            children = node_value['children']

            for child in children:

                self.networkx_graph.add_edge(child, parent_name)

        # Third Add the Parent Connections

        for node_key, node_value in self.network.items():

            parents = node_value['parents']
            child = node_value['name']

            for parent in parents:

                self.networkx_graph.add_edge(child, parent)

        return self.networkx_graph
