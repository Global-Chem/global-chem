#!/usr/bin/env python3
#
# GlobalChemExtensions - Development Operations
#
# ----------------------------------------------------

class ExtensionsError(Exception):

    __version_error_parser__ = "0.0.1"
    __allow_update__ = False

    '''
    
    Raise an Extension Error if something is wrong. 
    
    '''
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class QuantumChemistry(object):

    def __init__(self):

        self.name = 'quantum_chemistry'
        self.check_psi4_import()

    def check_psi4_import(self):

        '''

        Check to see if the psi4 component is installed

        '''

        try:
            import psi4
        except:
            print ("A Conda prequisite to use this package Psi4 not Installed")
            raise ImportError

    @staticmethod
    def get_zmatrix_store():

        '''
        
        Returns
            global_chem_zmatrix_store (Dict): Dictionary object of molecules to their respective zmatrix 

        '''

        from global_chem_extensions.quantum_chemistry.applications.zmatrix_store import ZMatrixStore

        zmatrix_store = ZMatrixStore()

        return zmatrix_store

    @staticmethod
    def get_psi4_parser():

        from global_chem_extensions.quantum_chemistry.applications.psi4_parser import Psi4Parser

        return Psi4Parser
