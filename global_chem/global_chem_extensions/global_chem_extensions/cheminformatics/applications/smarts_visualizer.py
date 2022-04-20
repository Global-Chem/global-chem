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

        self.smarts_pattern = self.prepare_smarts_for_url(smarts_pattern)


        self.base_url = "https://smarts.plus/smartsview/download_rest?"
        self.url = self.base_url + f'smarts={self.smarts_pattern}'

    def prepare_smarts_for_url(self, smarts_pattern):

        '''

        Replace some characters with URL accepted characters

        '''

        return smarts_pattern.replace("%", "%25").\
            replace("&", "%26").replace("+", "%2B").\
            replace("#", "%23").replace(";", "%3B")

    def get_image(self):

        res = requests.get(self.url)
        return Image(res.content)