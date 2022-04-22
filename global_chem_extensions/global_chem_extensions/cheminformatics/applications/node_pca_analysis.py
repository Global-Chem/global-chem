#!/usr/bin/env python3
#
# GlobalChemExtensions - Sunburster
#
# -------------------------------------

import sys

# Scientific Imports
# ------------------
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
colormaps = { 0: '#e6194b', 1: '#3cb44b',  2: '#ffe119', 3: '#4363d8', 4: '#f58231',  5: '#911eb4'}

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
                 save_file = False,
                 return_mol_ids = False,
        ):

        self.smiles_list = smiles_list
        self.morgan_radius = morgan_radius
        self.bit_representation = bit_representation
        self.number_of_clusters = number_of_clusters
        self.number_of_components = number_of_components
        self.random_state = random_state
        self.file_name = file_name
        self.save_file = save_file
        self.return_mol_ids = return_mol_ids

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
            x=chemicalspace[:,0],
            y=chemicalspace[:,2],
            img=[PCAAnalysis.mol2svg(m) for m in molecules_list],
            ids=[str(i) for i in range(0, len(self.smiles_list))],
            fill_color=kmeanc,
        )

        source = ColumnDataSource(kmean_data)
        plot = figure(plot_width=1000, plot_height=1000, tooltips=TOOLTIPS, title='Compounds')
        plot.circle('x', 'y',color='fill_color', size=10, fill_alpha=0.2,source=source)

        plot = gridplot([
            [plot]
        ])

        if self.save_file:

            output_file(filename=self.file_name, title="Static HTML file")
            save(plot)

        else:
            output_notebook()
            show(plot)

        if self.return_mol_ids:
            return self.plot_mappings



