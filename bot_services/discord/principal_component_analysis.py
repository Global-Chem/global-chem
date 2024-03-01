#!/usr/bin/env python3
#
# GlobalChemExtensions - Sunburster
#
# -------------------------------------

import sys

# Scientific Imports
# ------------------
import click
import numpy as np
import pandas as pd
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
               ):

    self.smiles_list = smiles_list
    self.morgan_radius = morgan_radius
    self.bit_representation = bit_representation
    self.number_of_clusters = number_of_clusters
    self.number_of_components = number_of_components
    self.random_state = random_state
    self.file_name = file_name
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

    output_file(filename=(self.file_name + ".html"), title="Static HTML file")
    save(plot)

# Click Commands
# --------------

@click.command()
@click.option('--smiles_list', default=[], help='Molecules to Analyze')
@click.option('--morgan_radius', default=1, help='Morgan Radius')
@click.option('--bit_representation', default=512, help='Bit Representation for the Fingerprint')
@click.option('--number_of_clusters', default=5, help='Number of Clusters to Categorize Data')
@click.option('--number_of_components', default=0.95, help='Number of Components')
@click.option('--random_state', default=0, help='')
@click.option('--file_name', default='pca_analysis.html', help='name of the html file')
@click.option('--principal_component_x', default=0, help='Principal Component of X-Axis')
@click.option('--principal_component_y', default=1, help='Principal Component of Y-Axis')
@click.option('--x_axis_label', default='PC1', help='X-Axis Title')
@click.option('--y_axis_label', default='PC2', help='Y-Axis Title')
@click.option('--plot_width', default=1000, help='Plot Width')
@click.option('--plot_height', default=1000, help='Plot Height')
@click.option('--title', default='Principal Component Analysis on SMILES', help='Title of Graph')

# Pipeline
# --------

def controller(
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
  plot_width,
  plot_height,
  title,
):

  '''

    Perform a pca analysis on a node within globalchem, can be extended to lists outside of the dedicated SMILES.

    Arguments:

        smiles_list (List): list of SMILES that the user wants to cluster
        morgan_radius (Int): Morgan Radius of the chemical environment
        bit_representation (Int): Length of the bit representation
        number_of_clusters (Int): Number of clusters the user would like to do
        number_of_components (Int): How many PCA vectors to analyze
        random_state (Int):
        file_name (String): file name the user would like to input
        save_file (Bool): Whether the user wants to display the plot or save it.
        return_mol_ids (Bool): Return the molecule IDS for the user to mine.
        principal_component_x (Int): Principal Component X on the x-axis that could be PC1 or PC2 or another vector
        principal_component_y (Int): Principal Component Y on the y-axis that could be PC1 or PC2 or another vector.
        x_axis_label (String): X Axis Label of the Plot
        y_axis_label (String): Y Axis Label of the Plot
        plot_width (Int): Width of the plot with default being 1000
        plot_height (Int): height of the plot
        title (String): Title of the plot

  '''

  _ = PCAAnalysis(
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
      plot_width,
      plot_height,
      title,
  )
if __name__ == '__main__':

    controller()



