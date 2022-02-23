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

        self.answer_layers = [ [], [] ]
        self.categories = {}
        self.dimensions = []

        self.total_values = {}

        self.labels = []
        self.parents = []
        self.values = []

        # Controller

        self._mine_global_chem()
        self._prepare_object()
        self.determine_matches()
        self.determine_possible_answers()
        self.create_answer_key_state()
        self.set_the_layers()
        self.radial_algorithm()
        self.sunburst()

    def _mine_global_chem(self):

        '''

        Mine the GlobalChem Data and get it ready for sunbursting

        '''

        gc = GlobalChem()
        nodes = gc.get_all_nodes()

        for node in nodes[1:]:

            node_initialized = node()

            node_name = node_initialized.name

            if node_name == 'common_regex_patterns': # Skip this
                continue

            node_smarts = node_initialized.get_smarts()

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

    def determine_matches(self):

        '''

        Determine the Substructure Matches using SMARTS with RDKit

        '''

        for row in self.compounds:

            smiles = row['smiles']
            molecule = Chem.MolFromSmiles(smiles)

            for index, chemical_list in enumerate(self.patterns):

                matches = []

                for name, pattern in chemical_list.items():

                    smarts_mol = Chem.MolFromSmarts(pattern)

                    try:
                        substructs = molecule.GetSubstructMatches(smarts_mol)
                        if substructs:
                            matches.append(name)
                    except:
                        pass

                row['%s:matches' % self.record_labels[index]] = matches

    def determine_possible_answers(self):

        '''

        Determine the possible answers

        '''

        possible_functional_groups = set()

        for row in self.compounds:
            for index, chemical_list in enumerate(self.patterns):
                functional_groups = row['%s:matches' % self.record_labels[index]]
                possible_functional_groups.update(functional_groups)

        if self.verbose:
            print ("Length of Possible Answers: %s" % len(possible_functional_groups))

        return possible_functional_groups

    def create_answer_key_state(self):

        '''

        Create the answer key state which is the layer, each list of lists is a layer

        '''

        self.answer_layers = [ [], [] ]

        for row in self.compounds:

            for index, chemical_list in enumerate(self.patterns):

                label = self.record_labels[index]
                functional_groups = row['%s:matches' % label]

                if len(functional_groups) > 0:
                    for j in functional_groups:
                        self.answer_layers[0].append(label)
                        self.answer_layers[1].append(j)


    def set_the_layers(self):

        '''

        Set the layers

        '''

        layers = [ self.answer_layers[0], self.answer_layers[1] ]
        questions = ['Record', 'Functional Group']

        for index, question in enumerate(questions):

            self.categories['label'] = question
            self.categories['values'] = layers[index]
            self.categories['categoryorder'] = 'category ascending'

            self.dimensions.append(self.categories)

    def radial_algorithm(self):

        '''

        Radial Algorithm, just go reread the plotly docs.

        '''


        for i in range(0, len(self.answer_layers[0])):

            paper = self.answer_layers[0][i]
            functional_group = self.answer_layers[1][i]
            key = paper + ' - ' + functional_group

            if key in self.total_values:
                self.total_values[key] += 1
            else:
                self.total_values[key] = 0

        # First Layer

        paper_values = {}

        for key , value in self.total_values.items():

            paper = key.split('-')[0]

            if paper not in paper_values:
                paper_values[paper] = 0
            else:
                paper_values[paper] += value


        first_layer_labels = ['Machine 1'] + self.record_labels
        first_layer_parents = [""] + ['Machine 1'] * len(self.record_labels)
        first_layer_values = [sum(paper_values.values())] + list(paper_values.values())

        # Second Layer

        second_layer_labels = [i for i in self.total_values.keys()]
        second_layer_parents = [i.split('-')[0].strip() for i in self.total_values.keys()]
        second_layer_values = [i for i in self.total_values.values()]

        self.labels = first_layer_labels + second_layer_labels
        self.parents = first_layer_parents + second_layer_parents
        self.values = first_layer_values + second_layer_values

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


