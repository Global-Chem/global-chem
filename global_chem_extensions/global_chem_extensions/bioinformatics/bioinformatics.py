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

class Bioinformatics(object):

    def __init__(self):

        self.name = 'bioinformatics'

    @staticmethod
    def initialize_globalchem_protein(
            pdb_id = None,
            pdb_path = None,
            peptide_sequence = None
    ):
        '''

        Arguments:
            pdb_id (String): pdb id unique 4 letter code
            pdb_path (String): path to the pdb file
            peptide_sequence (String): Peptide sequence
        '''

        from global_chem_extensions.bioinformatics.applications.protein import GlobalChemProtein

        global_chem_protein = GlobalChemProtein(
            pdb_file = pdb_path,
            fetch_pdb = pdb_id,
            peptide_sequence = peptide_sequence,
        )

        return global_chem_protein

    @staticmethod
    def initialize_globalchem_dna(
            dna_sequence,
            name = None,
    ):

        '''

        Arguments:
            dna_sequence (String): DNA sequence string
            name (String): Name of the DNA instance
        '''

        from global_chem_extensions.bioinformatics.applications.dna import GlobalChemDNA

        global_chem_dna = GlobalChemDNA(
            dna_sequence = dna_sequence,
            name = name
        )

        return global_chem_dna

    @staticmethod
    def initialize_globalchem_rna(
            rna_sequence,
            name = None,
    ):

        '''

        Arguments:
            rna_sequence (String): DNA sequence string
            name (String): Name of the DNA instance
        '''

        from global_chem_extensions.bioinformatics.applications.rna import GlobalChemRNA

        global_chem_rna = GlobalChemRNA(
            rna_sequence = rna_sequence,
            name = name
        )

        return global_chem_rna

    @staticmethod
    def initialize_globalchem_bacteria(
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

        from global_chem_extensions.bioinformatics.applications.bacteria import GlobalChemBacteria

        global_chem_bacteria = GlobalChemBacteria(
            name = name,
            genus_name = genus_name,
            species_name = species_name,
            family = family,
            cell_wall_composition = cell_wall_composition,
            gram_positive = gram_positive,
            gram_negative = gram_negative,
            inner_membrane = inner_membrane,
            outer_membrane = outer_membrane,
            shape = shape,
            nutrition_type = nutrition_type,
            respiratory = respiratory,
            environment = environment,
            antisera = antisera,
            universal_phylogenetic_tree_rna_sequence = universal_phylogenetic_tree_rna_sequence,
            taxonomic_hierarchy = taxonomic_hierarchy,
        )

        return global_chem_bacteria