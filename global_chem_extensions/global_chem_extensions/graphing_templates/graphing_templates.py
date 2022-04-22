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

class GraphingTemplates(object):

    def __init__(self):

        self.name = 'graphing_templates'

    @staticmethod
    def apply_plotly_template(
            figure,
            x_title='X-Axis',
            y_title = 'Y-Axis',
            height = 500,
            width = 1000
    ):
        '''

        Arguments:
            figure (Plotly Figure Object): plotly object you want beautified
            x_title (String): title of the x-axis
            y_title (String): title of the y-axis
            height (Int): height of the graph
            width (Int): width of the graph

        '''

        from  global_chem_extensions.graphing_templates.applications.plotly_template import PlotlyTemplate

        PlotlyTemplate(
            figure,
            x_title=x_title,
            y_title = y_title,
            height = height,
            width = width
        )
