#!/usr/bin/env python3
#
# GlobalChemExtensions - Dimorphite DL
#
# ------------------------------------

# Imports
# -------

from dimorphite_dl import DimorphiteDL

class DimorphiteAdapter(object):

    __version__ = '0.0.1'


    def __init__(self,
                 smiles,
                 min_ph,
                 max_ph,
                 max_variants,
                 label_states,
                pka_precision,
    ):

        '''

        Arguments:
            smiles (List): List of SMILES molecules

        '''

        self.min_ph = min_ph
        self.max_ph = max_ph
        self.pka_precision = pka_precision
        self.max_variants = max_variants
        self.label_states = label_states

        self.smiles = smiles
        self.output = {}

        self.dimorphite_dl = DimorphiteDL(
            min_ph=self.min_ph,
            max_ph=self.max_ph,
            max_variants=self.max_variants,
            label_states=self.label_states,
            pka_precision=self.pka_precision,
        )

    def run(self):

        for compound in self.smiles:

            self.output[compound] = self.dimorphite_dl.protonate(compound)

        return self.output