#!/usr/bin/env python3
#
# GlobalChemExtensions - Master Object
#
# -----------------------------------

class ExtensionsError(Exception):

    __version_error_parser__ = "0.0.1"
    __allow_update__ = False

    '''
    
    Raise an Extension Error if something is wrong. 
    
    '''
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class GlobalChemExtensions(object):

    def __init__(self):

        self.name = 'global_chem_extensions'

    @staticmethod
    def bioinformatics():

        '''

        Load the Development Operations Package

        '''

        from global_chem_extensions.bioinformatics.bioinformatics import Bioinformatics

        return Bioinformatics()

    @staticmethod
    def cheminformatics():

        '''

        Load the Cheminformatics Package

        '''

        from global_chem_extensions.cheminformatics.cheminformatics import ChemInformatics

        return ChemInformatics()

    @staticmethod
    def development_operations():

        '''

        Load the Development Operations Package

        '''

        from global_chem_extensions.development_operations.development_operations import DevelopmentOperations

        return DevelopmentOperations()

    @staticmethod
    def forcefields():

        '''

        Load the Development Operations Package

        '''

        from global_chem_extensions.forcefields.forcefields import ForceFields

        return ForceFields()

    @staticmethod
    def graphing_templates():

        '''

        Load the Graphing Templates Package

        '''

        from global_chem_extensions.graphing_templates.graphing_templates import GraphingTemplates

        return GraphingTemplates()

    @staticmethod
    def quantum_chemistry():

        '''

        Load the Quantum Chemistry Package

        '''

        from global_chem_extensions.quantum_chemistry.quantum_chemistry import QuantumChemistry

        return QuantumChemistry()