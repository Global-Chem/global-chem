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

class ForceFields(object):

    def __init__(self):

        self.name = 'force_fields'

    @staticmethod
    def initialize_globalchem_molecule(
            smiles,
            stream_file = None,
            frcmod_file = None,
    ):

        '''

        Arguments:
            smiles (String): A smiles string
            stream_file (String): stream file for CGenFF
            frcmod_file (String): FRCMOD file for GAFF2

        '''

        from global_chem_extensions.forcefields.applications.molecule import GlobalChemMolecule
        from global_chem_extensions.forcefields.applications.cgenff_molecule import CGenFFMolecule
        from global_chem_extensions.forcefields.applications.gaff2_molecule import GaFF2Molecule


        cgenff_molecule = None
        gaff2_molecule = None

        if stream_file:

            cgenff_molecule = CGenFFMolecule(stream_file=stream_file)

        if frcmod_file:

            gaff2_molecule = GaFF2Molecule(frcmod_file=frcmod_file)

        global_chem_molecule = GlobalChemMolecule(
            smiles=smiles,
            cgenff_molecule = cgenff_molecule,
            gaff2_molecule = gaff2_molecule
        )

        return global_chem_molecule

    @staticmethod
    def initialize_cgenff_molecule(stream_file):

        '''

        Arguments:

            stream_file (String): stream file from CGenFF

        '''

        from global_chem_extensions.forcefields.applications.cgenff_molecule import CGenFFMolecule

        cgenff_molecule = CGenFFMolecule(
            stream_file=stream_file
        )

        return cgenff_molecule

    @staticmethod
    def initialize_gaff2_molecule(frcmod_file):

        '''

        Arguments:

            frcmod_file (String): frcmod file for GaFF2

        '''

        from global_chem_extensions.forcefields.applications.gaff2_molecule import GaFF2Molecule

        gaff2_molecule = GaFF2Molecule(
            frcmod_file=frcmod_file
        )

        return gaff2_molecule

    @staticmethod
    def compute_cgenff_dissimilar_score(stream_file_1, stream_file_2, verbose=False):


        '''

        Arguments:
            stream_file_1 (String): file path of the first molecule
            stream_file_2 (String): file path of the second molecule
            verbose (Bool): If you want the full parameter output

        Returns:
            score (Float): similarity score

        '''

        from global_chem_extensions.forcefields.applications.cgenff_molecule import CGenFFMolecule
        from global_chem_extensions.forcefields.applications.dissimilarity_score import CGenFFDissimilarityScore

        molecule_1 = CGenFFMolecule(stream_file_1)
        molecule_2 = CGenFFMolecule(stream_file_2)

        dissimilar = CGenFFDissimilarityScore(
            molecule_1,
            molecule_2,
            verbose=verbose
        )

        score = dissimilar.compute_dissimilar_score_two_compounds()

        return score