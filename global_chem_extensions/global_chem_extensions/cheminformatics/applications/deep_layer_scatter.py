#!/usr/bin/env python3
#
# GlobalChemExtensions - DeepLayerScatter
#
# -------------------------------------------

# Plotly Imports & Kaleido Scopes
# -------------------------------

import plotly.express as px
import plotly.graph_objects as go
from kaleido.scopes.plotly import PlotlyScope

class DeepLayerScatter(object):

    '''

    Sunburster algorithm for analyzing a chemical list

    '''

    def __init__(self, deep_layer_network, save_file=False, verbose=False):

        self.deep_layer_network = deep_layer_network
        self.save_file = save_file

        self.verbose = verbose

        # Controller

        self.max_layer, self.node_counts = self.determine_max_layer_node_count()

    def determine_max_layer_node_count(self):

        '''

        Determine the Max Layer

        '''

        layer = 1
        node_counts = []

        for node_key, node_value in self.deep_layer_network.items():

            if int(node_value['layer']) > layer:
                layer = int(node_value['layer'])

        for i in range(1, layer + 1):

            counts = 0

            for node_key, node_value in self.deep_layer_network.items():

                if int(node_value['layer']) == i:
                    counts += 1

            node_counts.append(counts)

        return layer, node_counts

    def scatter(self, height = 800 , width = 1700):

        '''

        Scatter the Data

        Arguments:
            height (Int): Height of the graph
            width (Int): width of the graph

        TODO: Fix the values the weight later but get the essential visualization down.

        '''

        layer_object = []
        first_layer = [] # Establish the first layer to the root
        second_layer = [] # Establish the second layer to get the deep functionality going
        symmetric = False

        for i in range(1, self.max_layer + 1):

            nodes = [node_value for node_key, node_value in self.deep_layer_network.items() if int(node_value['layer']) == i]

            from functools import reduce

            if i == 1:
                first_layer = [ nodes[0]['name'] ] * reduce((lambda x, y: x * y), self.node_counts)
                layer_object.append(first_layer)

                values = nodes[0]['children']

                for value in values:
                    second_layer += [ value ] * self.node_counts[1]

                layer_object.append(second_layer)

            elif i ==2:
                pass
            else:

                layer = []

            for node in nodes:

                if not symmetric:
                    layer += [ node['name'] ]

                else:
                    layer += [ node['name'] ]  * self.node_counts[i - 1]

            if not symmetric:

                layer = layer  * self.node_counts[i - 1]

            layer_object.append(layer)

        layers_labels = [ str(i) for i in range(1, self.max_layer + 1 ) ]

        dimensions = []
        colors = []

        for index, layer in enumerate(layers_labels):

            categories = {}

            categories['label'] = layer
            categories['values'] = layer_object[index]

            if index == 1:
                colors += [ 'purple' ] * self.node_counts[0]

            elif index == 2:
                colors += [ 'lightseagreen' ] * self.node_counts[1]

            else:
                colors += [ 'gold' ] * self.node_counts[index - 1]


            categories['categoryorder'] = 'category ascending'

        dimensions.append(categories)

        fig = go.Figure(go.Parcats(
            dimensions=dimensions,
            arrangement='freeform',
            line = dict(
                color = colors,
            ),
            hoveron='color',
            hoverinfo='count+probability',
            labelfont={'size': 20, 'family': 'Times'},
            tickfont={'size': 15, 'family': 'Times'},
        ))

        fig.update_layout(
            title_text="Deep Layer Scatter",
            title_font=dict(size=12, family='Arial'),
            template='simple_white',
            xaxis_tickformat = 'i',
            bargap=0.4, # gap between bars of adjacent location coordinates,
            height=height,
            width=width
        )

        if self.save_file:

            ''' 
            
            Use this for large datasets
            
            '''

            scope = PlotlyScope(
                plotlyjs="https://cdn.plot.ly/plotly-latest.min.js",
            )
            with open("figure.png", "wb") as f:
                f.write(scope.transform(fig, format="png"))

        else:
            fig.show()




