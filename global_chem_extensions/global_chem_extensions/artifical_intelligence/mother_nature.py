#!/usr/bin/env python3
#
# GlobalChemExtensions - Artifical Intelligence
#
# --------------------------------------------

class ExtensionsError(Exception):

    __version_error_parser__ = "0.0.1"
    __allow_update__ = False

    '''
    
    Raise an Extension Error if something is wrong. 
    
    '''
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class MotherNature(object):

    def __init__(self):

        self.name = 'mother_nature'

        try:
            import tensorflow
        except:
            print ("A Conda prequisite to use this package Psi4 not Installed")
            raise ImportError


    @staticmethod
    def broca_area(
    ):

        from global_chem_extensions.mother_nature.applications.broca_area import BrocaArea

        return BrocaArea