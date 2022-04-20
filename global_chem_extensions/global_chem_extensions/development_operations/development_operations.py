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

class DevelopmentOperations(object):

    def __init__(self):

        self.name = 'development_operations'

    @staticmethod
    def check_status_on_open_source_databases():

        '''

        Check the Status on Databases

        '''

        from applications.database_monitor import DatabaseMonitor

        database_monitor = DatabaseMonitor()
        database_monitor.heartbeat()
