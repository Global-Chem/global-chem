#!/usr/bin/env python3
#
# GlobalChemExtensions - GAFF2 Molecule

# --------------------------------------

# Imports
# -------

import textwrap
import pandas as pd

class GaFF2Molecule(object):


    __version__ = '0.0.1'

    def __init__(self, frcmod_file, debugger=False):

        '''

        Arguments:
            gaff2_frcmod_file (String): File path of the gaff2 file

        '''

        # Debugging

        self.debugger = debugger

        # Conversion Variables

        self.counter = 1

        # Molecule Configurations

        self.atoms = []
        self.bond_connectivity = []
        self.improper_atoms = []
        self.residue_name = 'residue'
        self.molecule_charge = '0.000'

        # Parameters

        self.bond_parameters = []
        self.angle_parameters = []
        self.dihedral_parameters = []
        self.improper_dihedral_parameters = []
        self.nonbonded_parameters = []

        # Top Scores

        self.top_score_bond_parameter_penalty = 0
        self.top_score_charge_parameter_penalty = 0

        self.total_charge_penalty_score = 0

        self.total_bond_penalty_score = 0
        self.total_angle_penalty_score = 0
        self.total_dihedral_penalty_score = 0
        self.total_improper_penalty_score = 0

        self.frcmod_file = open(frcmod_file, 'r').read()
        self.lines = self.frcmod_file.split('\n')

        # Mine GAFF2 Data

        self.bonded_dataframe = pd.DataFrame()
        self.nonbonded_dataframe = pd.DataFrame()

        self.new_bonded_dataframe = pd.DataFrame()
        self.new_nonbonded_dataframe = pd.DataFrame()

        self.read_molecule_data()
        self.create_gaff2_dataframes()

    def read_molecule_data(self):

        '''

        Read the Molecule Data

        '''

        self.bond_parameters = list(filter(len, self.frcmod_file.split('BOND')[1].split('ANGLE')[0].split('\n')))
        self.angle_parameters = list(filter(len, self.frcmod_file.split('ANGLE')[1].split('DIHE')[0].split('\n')))
        self.dihedral_parameters = list(filter(len, self.frcmod_file.split('DIHE')[1].split('IMPROPER')[0].split('\n')))
        self.improper_dihedral_parameters = list(filter(len, self.frcmod_file.split('IMPROPER')[1].split('NONBON')[0].split('\n')))
        self.nonbonded_parameters = list(filter(len, self.frcmod_file.split('NONBON')[1].split('\n')))

        for row in self.nonbonded_parameters:
            if 'penalty score=' in row:
                self.total_charge_penalty_score += float(row.split('penalty score=')[1].strip().split()[0].strip(')'))

        for row in self.bond_parameters:
            if 'penalty score=' in row:
                self.total_bond_penalty_score += float(row.split('penalty score=')[1].strip().split()[0].strip(')'))

        for row in self.angle_parameters:
            if 'penalty score=' in row:
                self.total_angle_penalty_score += float(row.split('penalty score=')[1].strip().split()[0].strip(')'))

        for row in self.dihedral_parameters:
            if 'penalty score=' in row:
                self.total_dihedral_penalty_score += float(row.split('penalty score=')[1].strip().split()[0].strip(')'))

        for row in self.improper_dihedral_parameters:
            if 'penalty score=' in row:
                self.total_improper_penalty_score += float(row.split('penalty score=')[1].strip().split()[0].strip(')'))

    def create_gaff2_dataframes(self):


        '''

        Create the GAFF2 Dataframe:

            BOND
            atom type | force constant | bond length

            ANGLE
            atom type | force constant | angle degree

            DIHE
            atom type | idivf1 | force constant | phase | periodicity

            IMPROPER
            atom type | force constant | phase | periodicity

            NONBON
            atom type | rmin half | epsilon

        '''

        bond_parameter_types = []
        angle_parameter_types = []
        dihedral_parameter_types = []
        improper_parameter_types = []

        bonded_force_constants = []
        bonded_lengths = []

        angle_force_constants = []
        angle_degrees = []

        dihedral_fluctuations = []
        dihedral_force_constants = []
        dihedral_degrees = []
        dihedral_multipliticies = []

        improper_force_constants = []
        improper_degrees = []
        improper_multipliticies = []

        rmin_values = []
        epsilon_values = []
        nonbonded_parameter_types = []

        for row in self.nonbonded_parameters:

            atom_type = row.split()[0]
            rmin = row.split()[1]
            epsilon = row.split()[2]

            rmin_values.append(rmin)
            epsilon_values.append(epsilon)
            nonbonded_parameter_types.append(atom_type)

        for row in self.bond_parameters:

            row = row.lower()

            bond_type = (row.split('-')[0] + '-' + row.split('-')[1].split()[0]).replace(" ", "")

            bond_length = row.split('-')[1].split()[1]
            bond_force_constant =row.split('-')[1].split()[2]

            bond_parameter_types.append(bond_type)
            bonded_lengths.append(bond_length)
            bonded_force_constants.append(bond_force_constant)

        if self.debugger:
            print (bond_parameter_types)
            print (bonded_lengths)
            print (bonded_force_constants)

        for row in self.angle_parameters:

            row = row.lower()

            angle_type = (row.split('-')[0] + '-' + row.split('-')[1].split()[0] + '-' + row.split('-')[2].split()[0]).replace(" ", "")

            root = 'same'
            if 'using' in row:
                root = 'using'

            angle_degree = row.split(root)[0].split()[-1]
            angle_force_constant = row.split(root)[0].split()[-2]


            angle_degrees.append(angle_degree)
            angle_parameter_types.append(angle_type)
            angle_force_constants.append(angle_force_constant)

        if self.debugger:

            print (angle_degrees)
            print (angle_parameter_types)
            print (angle_force_constants)

        for row in self.dihedral_parameters:

            row = row.lower()

            dihedral_type = (row.split('-')[0] + '-' + row.split('-')[1].split()[0] + '-' + row.split('-')[2].split()[0] + '-' + row.split('-')[3].split()[0]).replace(" ", "")

            root = 'same'
            if 'using' in row:
                root = 'using'

            dihedral_fluctuation = row.split(root)[0].split()[-4]
            dihedral_force_constant = row.split(root)[0].split()[-3]
            dihedral_degree = row.split(root)[0].split()[-2]
            dihedral_multiplicity = row.split(root)[0].split()[-1]

            dihedral_parameter_types.append(dihedral_type)
            dihedral_fluctuations.append(dihedral_fluctuation)
            dihedral_force_constants.append(dihedral_force_constant)
            dihedral_degrees.append(dihedral_degree)
            dihedral_multipliticies.append(dihedral_multiplicity)

        if self.debugger:

            print (dihedral_parameter_types)
            print (dihedral_fluctuations)
            print (dihedral_force_constants)
            print (dihedral_degrees)
            print (dihedral_multipliticies)

        for row in self.improper_dihedral_parameters:

            row = row.lower()

            improper_type = '*' + (row.split('-')[0] + '-' + row.split('-')[1].split()[0] + '-' + row.split('-')[2].split()[0] + '-' + row.split('-')[3].split()[0]).replace(" ", "")

            root = 'same'
            if 'using' in row:
                root = 'using'

            improper_degree = row.split(root)[0].split()[-2]
            improper_multiplicity = row.split(root)[0].split()[-1]
            improper_force_constant = row.split(root)[0].split()[-3]

            improper_parameter_types.append(improper_type)
            improper_degrees.append(improper_degree)
            improper_multipliticies.append(improper_multiplicity)
            improper_force_constants.append(improper_force_constant)

        if self.debugger:

            print (improper_parameter_types)
            print (improper_degrees)
            print (improper_multipliticies)
            print (improper_force_constants)

        # atom type | idivf1 | force constant | phase | periodicity

        self.bonded_dataframe['atom_types'] = bond_parameter_types + angle_parameter_types + dihedral_parameter_types + improper_parameter_types
        self.bonded_dataframe['fluctuations'] = (len(bonded_lengths + angle_degrees) * [0]) + dihedral_fluctuations + (len(improper_degrees) * [0])
        self.bonded_dataframe['force_constants'] = bonded_force_constants + angle_force_constants + dihedral_force_constants + improper_force_constants
        self.bonded_dataframe['harmonic_value'] = bonded_lengths + angle_degrees + dihedral_degrees + improper_degrees
        self.bonded_dataframe['multiplicities'] = (len(bonded_lengths + angle_degrees) * [0]) + dihedral_multipliticies + improper_multipliticies

        self.nonbonded_dataframe['atom_types'] = nonbonded_parameter_types
        self.nonbonded_dataframe['rmin'] = rmin_values
        self.nonbonded_dataframe['epsilon'] = epsilon_values

        self.new_bonded_dataframe = self.bonded_dataframe.copy()
        self.new_nonbonded_dataframe = self.nonbonded_dataframe.copy()

    def update_bonded_dataframe(
            self,
            key,
            new_value,
            force_constant=True,
            harmonic_value=False,
            multiplicity=False,
            fluctuation=False,
    ):

        '''

        Update a value in the forcefield dataframes

        Arguments:

            key (String): Bonded key you would like to change
            new_value (String): new value to update
            replace_old_frame (Bool): Whether to replace the old frame
            force_constant (Bool): Whether the user wants to update the force constant
            harmonic_value (Bool): Whether the user wants to update the harmonic value
            multiplicity (Bool): Whether the user wants to update the mulicity of the rotation

        '''

        for i, row in self.new_bonded_dataframe.iterrows():

            if row['atom_types'] == key:

                if force_constant:

                    self.new_bonded_dataframe.at[i, 'force_constants' ] = new_value

                if harmonic_value:

                    self.new_bonded_dataframe.at[i, 'harmonic_value' ] = new_value

                if multiplicity:

                    self.new_bonded_dataframe.at[i, 'multiplicities' ] = new_value

                if fluctuation:
                    self.new_bonded_dataframe.at[i, 'fluctuations'] = new_value

    def update_nonbonded_dataframe(self, key, new_value, rmin=False, epsilon=False):

        '''

        Update a value in the forcefield dataframes

        Arguments:

            key (String): Bonded key you would like to change
            new_value (String): new value to update
            replace_old_frame (Bool): Whether to replace the old frame
            rmin (Bool): Lennard-Jones Parameter for the minimum equilbrium distance
            epsilon (Bool): Well-Depth of the energy potential for the R min.

        '''

        for i, row in self.new_nonbonded_dataframe.iterrows():

            if row['atom_types'] == key:

                if rmin:

                    self.new_bonded_dataframe.at[i, 'rmin' ] = new_value

                if epsilon:
                    self.new_bonded_dataframe.at[i, 'epsilon'] = new_value

    def write_frcmod_file(self, file_name = 'gaff2_modified.frcmod'):

        '''

        Write the GAFF2 FRCMOD File

        Arguments:
            file_name (String): File name of the gaff2 frcmod file

        '''

        out_file = open(file_name, 'w')

        gaff2_header = self._get_gaff2_header()

        # Header

        out_file.write(gaff2_header)

        # Bonds

        out_file.write("BONDS")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 2:
                out_file.write(f'{bond_parameter[0]}-{bond_parameter[1]}\t{row["force_constants"]}\t{row["harmonic_value"]}\n' )

        # Angles

        out_file.write('\n')
        out_file.write("ANGLE")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 3:
                out_file.write(f'{bond_parameter[0]}-{bond_parameter[1]}-{bond_parameter[2]}\t{row["harmonic_value"]}\t{row["force_constants"]}\n' )


        # DIHEDRALS

        out_file.write('\n')
        out_file.write("DIHE")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 4 and '*' not in row['atom_types']:
                out_file.write(f'{bond_parameter[0]}-{bond_parameter[1]}-{bond_parameter[2]}-{bond_parameter[3]}\t{row["fluctuations"]}\t{row["force_constants"]}\t{row["harmonic_value"]}\t{row["multiplicities"]}\n' )

        # IMPROPER DIHEDRALS

        out_file.write('\n')
        out_file.write("IMPROPER")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 4 and '*' in row['atom_types']:
                out_file.write(f'{bond_parameter[0].replace("*", "")}-{bond_parameter[1]}-{bond_parameter[2]}-{bond_parameter[3]}-{row["force_constants"]}\t{row["harmonic_value"]}\t{row["multiplicities"]} \n' )

        out_file.write('\n')
        out_file.write('\n')

        # Nonbond

        out_file.write("NONBON")
        out_file.write('\n')

        for i, row in self.nonbonded_dataframe.iterrows():

            nonbond_parameter = row['atom_types'].split('-')

            out_file.write(f'{nonbond_parameter[0]}-{nonbond_parameter[1]}\t{row["rmin"]}\t{row["epsilon"]}\n' )

    def _get_gaff2_header(self):

        '''

        Return the GAFF2 Header

        '''

        gaff2_header = '''
            Remark line goes here
            MASS
            '''

        return textwrap.dedent(gaff2_header)

    def __str__(self):

        return str(self.new_nonbonded_dataframe) + str(self.new_bonded_dataframe)