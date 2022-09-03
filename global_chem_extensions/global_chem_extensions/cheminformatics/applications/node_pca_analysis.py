#!/usr/bin/env python3
#
# GlobalChemExtensions - Sunburster
#
# -------------------------------------

import sys

# Scientific Imports
# ------------------
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

# RDkit Imports
# -------------

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import DataStructs
rdDepictor.SetPreferCoordGen(True)


# Graphing Imports
# ----------------

from bokeh.plotting import ColumnDataSource, figure, output_notebook, output_file, show, save
from bokeh.layouts import gridplot

# Global Configs
# --------------

TOOLTIPS = """<div>\nMolID: @ids<br>\n@img{safe}\n</div>\n"""
colormaps = {
    0: '#e6194b',
    1: '#3cb44b',
    2: '#ffe119',
    3: '#4363d8',
    4: '#f58231',
    5: '#911eb4',
    6: '#46f0f0', 
    7: '#f032e6', 
    8: '#bcf60c',
    9: '#fabebe', 
    10: '#008080',
    11: '#e6beff',
    12: '#9a6324', 
    13: '#fffac8', 
    14: '#800000',
    15: '#aaffc3',
    16: '#808000', 
    17: '#ffd8b1', 
    18: '#000075', 
    19: '#808080',
    20: '#ffffff', 
    21: '#000000'
}

class PCAAnalysis(object):

    __version__ = '0.0.1'

    def __init__(self,
                 smiles_list,
                 morgan_radius,
                 bit_representation,
                 number_of_clusters,
                 number_of_components,
                 random_state,
                 file_name,
                 principal_component_x,
                 principal_component_y,
                 x_axis_label,
                 y_axis_label,
                 plot_width ,
                 plot_height,
                 title,
                 save_file = False,
                 return_mol_ids = False,
                 save_principal_components = False,
        ):

        self.smiles_list = smiles_list
        self.morgan_radius = morgan_radius
        self.bit_representation = bit_representation
        self.number_of_clusters = number_of_clusters
        self.number_of_components = number_of_components
        self.random_state = random_state
        self.file_name = file_name
        self.save_file = save_file
        self.save_principal_components = save_principal_components
        self.return_mol_ids = return_mol_ids
        self.principal_component_x = principal_component_x
        self.principal_component_y = principal_component_y
        self.x_axis_label = x_axis_label
        self.y_axis_label = y_axis_label
        self.plot_width = plot_width
        self.plot_height = plot_height
        self.title = title

        self.plot_mappings = {}

    @staticmethod
    def mol2svg(mol):

        '''

        Convert the file from mol2 to svg

        Arguments:
            mol (RDKit): RDkit Mol Object

        Returns:
            d2d (DrawObject): object draw of the molecule.

        '''

        AllChem.Compute2DCoords(mol)
        d2d = rdMolDraw2D.MolDraw2DSVG(200,100)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def mol2fparr(self, mol):

        '''

        Convert the file from mol2 to fingerprint vector

        Arguments:
            mol (RDKit): RDkit Mol Object

        Returns:
            arr (List): Finger print ID

        '''

        arr = np.zeros((0,))
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, self.morgan_radius, nBits=self.bit_representation)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    def conduct_analysis(self):

        '''

        Conduct the PCA Analysis - comes from angels workflows for SMILES lists.

        '''

        molecules_list = [Chem.MolFromSmiles(i) for i in self.smiles_list]
        fingerprints_list = np.array([self.mol2fparr(m) for m in molecules_list])

        pca = PCA(n_components=self.number_of_components)
        chemicalspace = pca.fit_transform(fingerprints_list)
        kmean = KMeans(n_clusters=self.number_of_clusters, random_state=self.random_state)
        kmean.fit(fingerprints_list)
        kmeanc = [colormaps[i] for i in kmean.labels_]

        for i in range(0, len(self.smiles_list)):

            self.plot_mappings[str(i)] = self.smiles_list[i]

        kmean_data = dict(

            PCX = chemicalspace[:,self.principal_component_x],
            PCY = chemicalspace[:,self.principal_component_y],

            img=[PCAAnalysis.mol2svg(m) for m in molecules_list],
            ids=[str(i) for i in range(0, len(self.smiles_list))],
            fill_color=kmeanc,
        )

        ds= {
            'cols': [
                'ids','PCX','PCY','fill_color'
            ]
        }

        principal_component_df =pd.DataFrame(kmean_data)
        principal_component_df = principal_component_df[ds['cols']]
        principal_component_df['smiles'] = self.smiles_list

        source = ColumnDataSource(kmean_data)

        plot = figure(
            plot_width = self.plot_width,
            plot_height = self.plot_height,
            tooltips = TOOLTIPS,
            title = self.title,
            x_axis_label = self.x_axis_label,
            y_axis_label = self.y_axis_label,
        )
        plot.circle('PCX','PCY', color='fill_color', size=10, fill_alpha=0.2, source=source)

        plot = gridplot([
            [plot]
        ])

        if self.save_file:

            output_file(filename=(self.file_name + ".html"), title="Static HTML file")
            save(plot)

        else:
            output_notebook()
            show(plot)
            
        if self.save_principal_components:
            principal_component_df.to_csv((self.file_name + ".tsv"), sep="\t", index=False)

        if self.return_mol_ids:
            return self.plot_mappings



