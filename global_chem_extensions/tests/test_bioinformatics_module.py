# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

# imports
# -------

from global_chem_extensions import GlobalChemExtensions


def test_load_bioinformatics_module():

    '''

    Test the Loading of the Bioinformatics Module

    '''

    gce = GlobalChemExtensions()
    gce.bioinformatics()