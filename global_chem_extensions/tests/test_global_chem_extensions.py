# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

# imports
# -------

from global_chem_extensions import GlobalChemExtensions

def test_initialize_class():

    '''

    Test the intialize of the class, with parameters if it extends to that

    '''

    GlobalChemExtensions()

def test_load_bioinformatics_module():

    '''

    Test the Loading of the Bioinformatics Module

    '''

    gce = GlobalChemExtensions()
    gce.bioinformatics()

def test_load_cheminformatics_module():

    '''

    Test the Loading of the Cheminformatics Module

    '''

    gce = GlobalChemExtensions()
    gce.cheminformatics()

def test_load_development_operations_module():

    '''

    Test the Loading of the Development Operations Module

    '''

    gce = GlobalChemExtensions()
    gce.development_operations()

def test_load_quantum_chemistry_module():

    '''

    Test the Loading of the Quantum Chemistry Module

    '''

    gce = GlobalChemExtensions()
    # gce.quantum_chemistry()

def test_load_forcefields_module():

    '''

    Test the Loading of the ForceFields Module

    '''

    gce = GlobalChemExtensions()
    gce.forcefields()