#!/usr/bin/env python3
#
# GlobalChemExtensions - Chemical Diversity

# -----------------------------------------

# Imports
# -------

import pandas as pd
from rdkit import Chem
import plotly.graph_objects as go
from global_chem import GlobalChem
from rdkit.Chem import DataStructs

if __name__ == '__main__':

  gc = GlobalChem()
  gc.build_global_chem_network()
  node_items = gc.get_node_smiles('pihkal')



