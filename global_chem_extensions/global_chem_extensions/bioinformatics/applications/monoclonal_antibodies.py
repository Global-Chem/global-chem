#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Monoclonal Antibodies

# -------------------------------------------------------

class GlobalChemMonoclonalAntibody(object):

    __COMMON_ANTIBODIES__ = [
        'adalimumab_humira',
        'basiliximab_simulect',
        'infliximab_avsola',
        'infliximab_inflextra',
        'infliximab_infliximab',
        'infliximab_remicade',
        'infliximab_renflexis',
        'palivizumab_synagis',
        'rituximab_riabni',
        'rituximab_rituxan',
        'rituximab_rituxan_hycela',
        'rituximab_ruxience',
        'rituximab_truxima',
        'trastuzumab_enhertu',
        'trastuzumab_herceptin',
        'trastuzumab_herzuma',
        'trastuzumab_kadcyla',
        'trastuzumab_kanjinti',
        'trastuzumab_ogivri',
        'trastuzumab_ontruzant'
    ]

    def __init__(self, name):

        self.name = name

    def get_common_antibodies(self):

        '''

        Returns the Common antibodies

        '''

        return self.__COMMON_ANTIBODIES__