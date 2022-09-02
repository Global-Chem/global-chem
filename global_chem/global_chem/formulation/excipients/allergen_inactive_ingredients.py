#!/usr/bin/env python3
#
# GlobalChem - Allergen Inactive Ingredients
# Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7122736/
# ------------------------------------------

class AllergenInactiveIngredients(object):

    def __init__(self):

        self.name = 'allergen_inactive_ingredients'

    @staticmethod
    def get_smiles():

        smiles = {
            'lactose': '',
            'corn starch': '',
            'povidone': '',
            'carboxymethylcellulose': '',
            'gelatin': '',
            'brilliant blue	dye': '',
            'sunset yellow fcf': '',
            'allura red	dye': '',
            'propylene glycol': '',
            'indigo carmine	dye': '',
            'mannitol': '',
            'sucrose': '',
            'sodium benzoate': '',
            'parabens': '',
            'aspartame': '',
            'erythrosine dye': '',
            'tartrazine dye': '',
            'saccharin': '',
            'poloxamer': '',
            'soybean oil': '',
            'benzyl alcohol': '',
            'vanilla': '',
            'castor oil': '',
            'cetyl alcohol': '',
            'sulfite': '',
            'polyethylene glycol castor oil': '',
            'peanut oil': '',
            'benzoic acid': '',
            'corn syrup': '',
            'sesame oil': '',
            'starch wheat': '',
            'casein': '',
            'banana essence': '',
            'milk': '',
            'glucosamine': '',
            'new coccine dye': '',
            'stearyl alcohol': '',

        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {

        }

        return smarts

    @staticmethod
    def get_bit_vector():

        bit_vector = {

        }

        return bit_vector