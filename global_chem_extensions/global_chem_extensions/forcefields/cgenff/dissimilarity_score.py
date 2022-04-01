#!/usr/bin/env python3
#
# GlobalChemExtensions - Dissimilarity Score

# ------------------------------------------

from .cgenff_molecule import CGenFFMolecule

class CGenFFDissimilarityScore(object):

    __version__ = '0.0.1'

    '''
    
    Here's the first draft of the score: 
    
    We take into account SMART pattern recognition between two functional groups. 
    
    '''

    def __init__(self, molecule_1, molecule_2, verbose=False):


        '''
        Arguments:
            molecule_1 (CGenFFMolecule Object): List of SMILES strings
            molecule_2 (CGenFFMolecule Object): List of SMILES strings
        '''

        self.molecule_1 = molecule_1
        self.molecule_2 = molecule_2
        self.verbose = verbose

    def _compute_bonded_term_similarity(self):

        '''

        Compute the Similarity between the two compounds bonded terms.

        '''

        bond_weight = 4

        penalty_rows_1 = self.molecule_1.bond_parameters
        penalty_rows_2 = self.molecule_2.bond_parameters

        bond_terms_penalty_row_1 = [ (i.split()[0],i.split()[1]) for i in penalty_rows_1 ]
        bond_terms_penalty_row_2 = [ (i.split()[0],i.split()[1]) for i in penalty_rows_2 ]

        bond_terms_identical = set( bond_terms_penalty_row_1 ) & set(bond_terms_penalty_row_2 )
        bond_terms_difference_1 = sum( ([x,y] for x,y in zip( bond_terms_penalty_row_1, bond_terms_penalty_row_2 ) if x != y ), [])

        total_terms = len(bond_terms_penalty_row_1) + len(bond_terms_penalty_row_2)
        score = total_terms - len(bond_terms_identical)

        if len(bond_terms_difference_1) == 0:

            dissimilar_score = 0

        else:

            dissimilar_score = (len(bond_terms_difference_1) - len(bond_terms_identical)) * bond_weight

        dissimilar_percentage = 100 * (1 - (len(bond_terms_identical) / score))

        if self.verbose:
            print ("Bond Dissimilar Percentage: %s " % dissimilar_percentage)

        return dissimilar_score

    def _compute_angle_term_similarity(self):

        '''
        Compute the Similarity between the two compounds angle terms.
        '''

        angle_weight = 3

        penalty_rows_1 = self.molecule_1.angle_parameters
        penalty_rows_2 = self.molecule_2.angle_parameters

        angle_terms_penalty_row_1 = [ (i.split()[0],i.split()[1], i.split()[2]) for i in penalty_rows_1 ]
        angle_terms_penalty_row_2 = [ (i.split()[0],i.split()[1], i.split()[2]) for i in penalty_rows_2 ]

        angle_terms_identical = set( angle_terms_penalty_row_1) & set(angle_terms_penalty_row_2 )
        angle_terms_difference = sum( ([x,y] for x,y in zip( angle_terms_penalty_row_1, angle_terms_penalty_row_2 ) if x != y ), [])


        total_terms = len(penalty_rows_1) + len(penalty_rows_2)
        score = total_terms - len(angle_terms_identical)
        dissimilar_percentage = 100 * (1 - (len(angle_terms_identical) / score))

        if self.verbose:
            print ("Angle Dissimilar Percentage: %s " % dissimilar_percentage)

        dissimilar_score = (len(angle_terms_difference) - len(angle_terms_identical)) * angle_weight


        return dissimilar_score

    def _compute_dihedral_term_similarity(self):

        '''

        Compute the Similarity between the two compounds dihedral terms.

        '''

        dihedral_weight = 2

        penalty_rows_1 = self.molecule_1.dihedral_parameters
        penalty_rows_2 = self.molecule_2.dihedral_parameters

        dihedral_terms_penalty_row_1 = [ (i.split()[0],i.split()[1], i.split()[2], i.split()[3]) for i in penalty_rows_1 ]
        dihedral_terms_penalty_row_2 = [ (i.split()[0],i.split()[1], i.split()[2], i.split()[3]) for i in penalty_rows_2 ]

        dihedral_terms_identical = set( dihedral_terms_penalty_row_1) & set(dihedral_terms_penalty_row_2 )
        dihedral_terms_difference = sum( ([x,y] for x,y in zip( dihedral_terms_penalty_row_1, dihedral_terms_penalty_row_2 ) if x != y ), [])

        total_terms = len(penalty_rows_1) + len(penalty_rows_2)

        score = total_terms - len(dihedral_terms_identical)
        dissimilar_score = (len(dihedral_terms_difference) - len(dihedral_terms_identical)) * dihedral_weight

        dissimilar_percentage = 100 * (1 - (len(dihedral_terms_identical) / score))

        if self.verbose:
            print ("Dihedral Dissimilar Percentage: %s " % dissimilar_percentage)

        return dissimilar_score

    def _compute_improper_dihedral_term_similarity(self):

        '''

        Compute the Similarity between the two compounds dihedral terms.

        '''

        improper_dihedral_weight = 2

        penalty_rows_1 = self.molecule_1.improper_dihedral_parameters
        penalty_rows_2 = self.molecule_2.improper_dihedral_parameters

        improper_dihedral_terms_penalty_row_1 = [ (i.split()[0],i.split()[1], i.split()[2], i.split()[3]) for i in penalty_rows_1 ]
        improper_dihedral_terms_penalty_row_2 = [ (i.split()[0],i.split()[1], i.split()[2], i.split()[3]) for i in penalty_rows_2 ]

        improper_dihedral_terms_identical = set( improper_dihedral_terms_penalty_row_1) & set(improper_dihedral_terms_penalty_row_2 )
        improper_dihedral_terms_difference = sum( ([x,y] for x,y in zip( improper_dihedral_terms_penalty_row_1, improper_dihedral_terms_penalty_row_2 ) if x != y ), [])

        dissimilar_score = (len(improper_dihedral_terms_difference) - len(improper_dihedral_terms_identical)) * improper_dihedral_weight

        total_terms = len(penalty_rows_1) + len(penalty_rows_2)
        score = total_terms - len(improper_dihedral_terms_identical)
        dissimilar_percentage = 100 * (1 - (len(improper_dihedral_terms_identical) / score))

        if self.verbose:
            print ("Improper Dihedral Dissimilar Percentage: %s " % dissimilar_percentage)

        return dissimilar_score

    def compute_dissimilar_score_two_compounds(self):

        '''

        Compute the dissimilarity between two compounds

        '''

        compound_1 = self.molecule_1
        compound_2 = self.molecule_2

        bond_dissimilar_score = abs(self._compute_bonded_term_similarity())
        angle_dissimilar_score = abs(self._compute_angle_term_similarity())
        dihedral_dissimilar_score = abs(self._compute_dihedral_term_similarity())
        improper_dihedral_dissimilar_score = abs(self._compute_improper_dihedral_term_similarity())

        dissimilar_score = (
                abs(bond_dissimilar_score) +
                abs(angle_dissimilar_score) +
                abs(dihedral_dissimilar_score) +
                abs(improper_dihedral_dissimilar_score)
        )

        if self.verbose:
            print ( "Bond Dissimilar Score: %s " % bond_dissimilar_score )
            print ( "Angle Dissimilar Score: %s " % angle_dissimilar_score)
            print ( "Dihedral Dissimilar Score: %s " % dihedral_dissimilar_score)
            print ( "Improper Dihedral Dissimilar Score: %s " % improper_dihedral_dissimilar_score)

        return dissimilar_score