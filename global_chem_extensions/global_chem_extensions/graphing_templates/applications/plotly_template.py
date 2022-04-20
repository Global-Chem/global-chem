#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Plotly Template

# --------------------------------------------------


class PlotlyTemplate(object):


    __version__ = '0.0.1'

    def __init__(
            self,
            figure,
            x_title= 'x-axis',
            y_title = 'y-axis',
            height = 500,
            width = 1000
    ):

        self.figure = figure
        self.x_title = x_title
        self.y_title = y_title
        self.height = height
        self.width = width

        self.apply_template()

    def apply_template(self):

        '''

        Apply the template

        '''

        self.figure.update_layout(legend=dict(itemsizing='constant'))
        self.figure.update_layout(legend=dict(
            orientation="v",
            yanchor="bottom",
            y=0.96,
            xanchor="right",
            x=1,
            font = dict(family = "Arial", size = 10),
            bordercolor="LightSteelBlue",
            borderwidth=2,
        ),
            legend_title = dict(font = dict(family = "Arial", size = 10))
        )

        # ----- Handle Layout -----

        self.figure.update_xaxes(
            ticks="outside",
            tickwidth=1,
            tickcolor='black',
            tickfont=dict(family='Arial', color='black', size=13),
            title_font=dict(size=13, family='Arial'),
            title_text=self.x_title,
            ticklen=15,
        )

        self.figure.update_yaxes(
            ticks="outside",
            tickwidth=1,
            tickcolor='black',
            title_text=self.y_title,
            tickfont=dict(family='Arial', color='black', size=13),
            title_font=dict(size=13, family='Arial'),
            ticklen=15,

        )

        self.figure.update_layout(
            xaxis_tickformat = 'i',
            bargap=0.2, # gap between bars of adjacent location coordinates,
            height=self.height,
            width=self.width,
            plot_bgcolor='rgba(0,0,0,0)'
        )