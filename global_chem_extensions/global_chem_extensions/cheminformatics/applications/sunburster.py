#!/usr/bin/env python3
#
# GlobalChemExtensions - Sunburster
#
# -------------------------------------

# RDKit Imports
# -------------

from rdkit import Chem

# GlobalChem Imports
# ------------------

from global_chem import GlobalChem

# Plotly Imports & Kaleido Scopes
# -------------------------------

import plotly.graph_objects as go
from kaleido.scopes.plotly import PlotlyScope

class Sunburster(object):

    '''

    Sunburster algorithm for analyzing a chemical list

    '''

    def __init__(self, smiles_list, save_file=False, verbose=False):

        self.smiles_list = smiles_list
        self.save_file = save_file

        self.patterns = []
        self.record_labels = []

        self.compounds = []
        self.verbose = verbose

        self.category_first_layer = []
        self.category_counts = {}

        self.total_values = {}

        self.labels = []
        self.parents = []
        self.values = []

        # Second Layer

        self.second_labels = []
        self.second_parents = []
        self.second_values = []

        # Controller

        self._mine_global_chem()
        self._prepare_object()
        self.determine_first_layer()
        self.determine_second_layer()
        self.set_the_layers()
        self.sunburst()

    def _mine_global_chem(self):

        '''

        Mine the GlobalChem Data and get it ready for sunbursting

        '''

        gc = GlobalChem()
        nodes = gc.get_all_nodes()

        for node in nodes[1:]:

            node_name = node.name

            if node_name == 'common_regex_patterns': # Skip this
                continue

            node_smarts = node.get_smarts()

            self.record_labels.append(node_name)
            self.patterns.append(node_smarts)

    def _prepare_object(self):

        '''

        Prepare the object of the smiles, just get it into the format.

        '''

        for smiles in self.smiles_list:

            row = {}
            row['smiles'] = smiles
            self.compounds.append(row)

    def determine_first_layer(self):

        '''

        Determine the Record Category for the First Layer of the Sunburst

        '''

        for row in self.compounds:

            smiles = row['smiles']
            molecule = Chem.MolFromSmiles(smiles)

            for index, chemical_list in enumerate(self.patterns):

                matches = []

                for name, pattern in chemical_list.items():

                    smarts_mol = Chem.MolFromSmarts(pattern)

                    try:
                        substructs = molecule.HasSubstructMatch(smarts_mol)
                        if substructs:
                            matches.append(name)
                    except:
                        pass

                longest_name = max(matches, key=str)
                shortest_name = min(matches, key=str)

                if longest_name not in self.category_first_layer:
                    self.category_first_layer.append(longest_name)

                if shortest_name not in self.category_first_layer:
                    self.category_first_layer.append(shortest_name)

                if longest_name not in self.category_counts:
                    self.category_counts[longest_name] = 1
                else:
                    self.category_counts[longest_name] += 1

                if shortest_name not in self.category_counts:
                    self.category_counts[shortest_name] = 1
                else:
                    self.category_counts[shortest_name] += 1

                row['%s:matches' % self.record_labels[index]] = matches

                row['matches'] = [ longest_name , shortest_name ]
                row['relation'] = longest_name + ' - ' + shortest_name

    def determine_second_layer(self):

        '''

        Determine the Substructure Matches using SMARTS with RDKit for the Category Second Layer

        '''

        for compound in self.compounds:

            if 'matches' in compound:

                parent = compound['matches'][0]
                label = compound['relation']

                self.second_parents.append(parent)

                if label not in self.second_labels:
                    self.second_labels.append(label)
                self.second_values.append(1)

        for i in range(0, len(self.second_labels)):

            parent = self.second_labels[i].split(' - ')[0].strip()

            if parent not in self.total_values:
                self.total_values[parent] = 0

            self.total_values[parent] += self.second_values[i]

    def set_the_layers(self):

        '''

        Set the layers

        '''

        category_first_layer_labels = ['GlobalChem'] + self.category_first_layer
        category_first_layer_parents = [""] + ["Globalchem"] * len(self.category_first_layer)
        category_first_layer_values = [sum(self.total_values.values())] + list(self.total_values.values())

        category_second_layer_labels = self.second_labels
        category_second_layer_parents = self.second_parents
        category_second_layer_values = self.second_values

        # Construct the Layers

        self.labels = category_first_layer_labels + category_second_layer_labels
        self.parents = category_first_layer_parents + category_second_layer_parents
        self.values = category_first_layer_values + category_second_layer_values

    def sunburst(self):

        '''

        Sunburst the Data

        '''

        fig =go.Figure()

        fig.add_trace(go.Sunburst(
            labels=self.labels,
            parents=self.parents,
            values=self.values,
            insidetextorientation='radial',
            #     marker={"colors": colors_charge},
            domain=dict(column=0),
            opacity=0.9,
        ))

        fig.update_layout(
            margin = dict(t=5, l=5, r=5, b=5),
            height=600,
            width=600,
            grid= dict(columns=1, rows=1)
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