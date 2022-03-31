#!/usr/bin/env python3
#
# GlobalChemExtensions - SMARTS Pattern Visualizer
#
# ------------------------------------------------


# Imports
# -------

import urllib
import requests

from rdkit import Chem
from IPython.display import Image

class SmartsVisualizer(object):

    __version__ = '0.0.1'


    def __init__(self,
                 smarts_pattern
        ):

        self.smarts_pattern = smarts_pattern

        self.base_url = "https://smarts.plus/smartsview/download_rest?"
        self.url = self.base_url + f"smarts={self.smarts_pattern}"

    def get_image(self):

        res = requests.get(self.url)
        return Image(res.content)