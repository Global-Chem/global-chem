#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Bacteria

# -----------------------------------------

class GlobalChemBacteria(object):

    __ALLOWED_CELL_WALL_COMPOSITIONS__ = [
        'thick_peptidoglycan',
        'lipopolysaccharide'
    ]

    __BACTERIA_SHAPES__ = {
        'bacillus': 'rod',
        'spirilla': 'spiral',
        'coccus': 'sphere',
        'vibrio': 'comma',
    }

    __NUTRITIONAL_SOURCES__ = {
        'anoxygenic_photosynthetic' : 'thiosulfate',
        'autotrophic': 'self',
        'chemoautotrophs': 'inorganic_oxidation',
        'green_sulfur': 'sulfide:hydyogen:iron',
        'oxygenic_photosynthetic': 'oxygenic_photosynthesis',
        'photoautotrophs': 'light',
        'purple_sulfur': 'thiosulphates',
        'purple_non_sulfur': 'hydrogen',
        'sulfur': 'hydrogen_sulfide:sulphate:thiosulphate',
        'nitrogen': 'nitrite:nitrate',
        'hydrogen': 'oxidize_hydrogen',
        'methanotrophs': 'methane',
        'iron': 'oxidize_ferrous_ions',
        'heterotophic': 'organic_compounds',
        'parasitic': 'living_organisms',
        'symbiotic': 'symbiotic_association',
        'saprophytic': 'decaying_organic_matter',
    }

    __RESPIRATION_METHODS__ = {
        'aerobic': ' molecular_oxygen',
        'obligate': 'only_aerobically',
        'facultative_anaerobes': 'anaerobically:anaerobic',
        'anaerobic': 'carbon_dioxide:sulphur:fumarate:iron_three',
    }

    __ENVIRONMENT_RESERVOIRS__ = {
        'endogenous': 'inside_body',
        'exogenous': 'outside_body',
    }

    __SEROLOGIC_SYSTEMS__ = [
        'carbohydrate',
        'protein'
    ]

    __TAXONOMIC_CLASSIFIERS__ = [
        'domain',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species'
    ]

    def __init__(
            self,
            name = 'staphylococcus',
            genus_name = '',
            species_name = '',
            family = 'chlorobiaceae',
            cell_wall_composition = 'thick_peptidoglycan',
            gram_positive = False,
            gram_negative = False,
            inner_membrane = None,
            outer_membrane = None,
            shape = None,
            nutrition_type = None,
            respiratory = None,
            environment = None,
            antisera = None,
            universal_phylogenetic_tree_rna_sequence = None,
            taxonomic_hierarchy = None,

    ):

        '''

        Arguments:
            name (String): Name of the Bacteria
            family (String): Family of the Bacteria
            genus_name (String): Genus Name of the Bacteria
            species_name (String): Species Name of the Bacteria where it is derived from
            cell_wall_composition (String): What the composition of the cell wall is made up of
            gram_positive (Bool): Typically means one plasma membrane and set as a Bool Flag for user use.
            gram_negative (Bool): Typically means two plasma membrane and set as a Bool Flag for user use.
            inner_membrane (String): Seperates the Cytoplasm from the periplasm
            outer_membrane (String): Seperates the Periplasm from the extracellular space
            shape (String): Classification based on shape and needs comparison to the allowed shape classifiers
            nutrition_type (String): the nutritional source of the bacteria to gain energy.
            environment (String): environment of the bacteria
            antisera (String): the Serologic environment TODO: What does this mean?
            universal_phylogenetic_tree_rna_sequence (Dict): RNA Sequence
            taxonomic_hierarchy (Dict): A hierarchy of the taxonomy

        References:
            1.) https://byjus.com/biology/bacteria/#classification
            2.) https://byjus.com/neet/autotrophic-bacteria/
            3.) https://byjus.com/neet/heterotrophic-bacteria/
            4.) https://byjus.com/biology/difference-between-aerobic-and-anaerobic-bacteria/
            5.) https://www.biologydiscussion.com/bacteria/respiration-in-bacteria-with-diagram/52066
            6.) http://www.columbia.edu/itc/hs/medical/pathophys/id/2009/introNotes.pdf

        '''


        # All Variables
        # --------------

        self.name = name
        self.family = family
        self.genus_name = genus_name
        self.species_name = species_name

        # Phenotype Classifiers
        # ---------------------


        # Based on Cell Wall and Gram-Negative and Gram-Positive

        self.cell_wall_composition = cell_wall_composition
        self.gram_positive = gram_positive
        self.gram_negative = gram_negative
        self.inner_membrane = inner_membrane
        self.outer_membrane = outer_membrane

        # Based on Shape

        if shape and shape not in list(self.__BACTERIA_SHAPES__.values()):
            print ("Shape not in the accepted bacteria shape classifiers")
            raise ValueError
        else:
            self.shape = shape

        # Based on Nutrition

        if nutrition_type and nutrition_type not in list(self.__NUTRITIONAL_SOURCES__.values()):
            print ("Nutrition not in the accepted bacteria shape classifiers")
            raise ValueError
        else:
            self.nutrition_type = nutrition_type

        # Based on Respiratory

        if respiratory and respiratory not in list(self.__RESPIRATION_METHODS__.values()):
            print ("Respiratory not in the accepted bacteria shape classifiers")
            raise ValueError
        else:
            self.nutrition_type = nutrition_type

        # Based on Environment

        if environment and environment not in list(self.__ENVIRONMENT_RESERVOIRS__.values()):
            print ("Environment not in the accepted bacteria shape classifiers")
            raise ValueError
        else:
            self.environment = environment

        # Based on Serologic

        if antisera and antisera not in list(self.__SEROLOGIC_SYSTEMS__):
            print ("Antisera not in the accepted bacteria shape classifiers")
            raise ValueError
        else:
            self.antisera = antisera

        # Based on Tree

        self.universal_phylogenetic_tree_rna_sequence = universal_phylogenetic_tree_rna_sequence

        # Based on Taxonomic Classifiers

        if taxonomic_hierarchy:
            for classifier in self.__TAXONOMIC_CLASSIFIERS__:
                if classifier not in taxonomic_hierarchy:
                    print ("Missing Classifier: %s" % classifier)
                    raise ValueError

        self.taxonomic_hierarchy = taxonomic_hierarchy

    def get_classifiers(self):

        '''

        Get all the classifiers for the bacteria object

        Returns:
            classifiers (Dict): Dictionary of the bacteria classifiers

        '''

        classifiers = {
            'allowed_cell_wall_compositions': self.__ALLOWED_CELL_WALL_COMPOSITIONS__,
            'bacteria_shapes': self.__BACTERIA_SHAPES__,
            'nutritional_sources': self.__NUTRITIONAL_SOURCES__,
            'respiratory_methods': self.__RESPIRATION_METHODS__,
            'environment_reservoirs': self.__ENVIRONMENT_RESERVOIRS__,
            'serologic_systems': self.__SEROLOGIC_SYSTEMS__,
            'taxonomic_classifiers': self.__TAXONOMIC_CLASSIFIERS__
        }

        return classifiers