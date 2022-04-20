#!/usr/bin/env python3
#
# GlobalChemExtensions - CGenFF Molecule

# --------------------------------------

# Imports
# -------

import textwrap
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import GetPeriodicTable

class CGenFFMolecule(object):


    __version__ = '0.0.1'

    def __init__(self, stream_file):

        '''

        Arguments:
            cgenff_stream_file (String): File path of the cgenff file

        '''

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

        # Top Scores

        self.top_score_bond_parameter_penalty = 0
        self.top_score_charge_parameter_penalty = 0

        self.total_charge_penalty_score = 0

        self.total_bond_penalty_score = 0
        self.total_angle_penalty_score = 0
        self.total_dihedral_penalty_score = 0
        self.total_improper_penalty_score = 0

        self.stream_file = open(stream_file, 'r').read()
        self.lines = self.stream_file.split('\n')

        # Mine CGenFF Data

        self.bonded_dataframe = pd.DataFrame()
        self.nonbonded_dataframe = pd.DataFrame()

        self.new_bonded_dataframe = pd.DataFrame()
        self.new_nonbonded_dataframe = pd.DataFrame()

        self.read_molecule_data()
        self.create_cgenff_dataframes()

    def read_molecule_data(self):

        '''

        Read the Molecule Data

        '''

        self.atoms = list(filter(lambda x: 'ATOM' in x, self.lines))
        self.bond_connectivity = list(filter(lambda x: 'BOND ' in x, self.lines))
        self.improper_atoms = list(filter(lambda x: 'IMPR ' in x, self.lines)) # Don't get confused with IMPROPERS add space

        self.bond_parameters = list(filter(len, self.stream_file.split('BONDS')[1].split('ANGLES')[0].split('\n')))
        self.angle_parameters = list(filter(len, self.stream_file.split('ANGLES')[1].split('DIHEDRALS')[0].split('\n')))
        self.dihedral_parameters = list(filter(len, self.stream_file.split('DIHEDRALS')[1].split('IMPROPERS')[0].split('\n')))
        self.improper_dihedral_parameters = list(filter(len, self.stream_file.split('IMPROPERS')[1].split('END')[0].split('\n')))

        penalty_row = list(filter(lambda x: 'param penalty' in x, self.lines))[0]

        self.top_score_bond_parameter_penalty = penalty_row.split('param penalty=')[1].split(';')[0].strip()
        self.top_score_charge_parameter_penalty = penalty_row.split('charge penalty=')[1].strip()
        self.residue_name = penalty_row.split('RESI')[1].strip().split()[0].strip()
        self.molecule_charge = penalty_row.split('RESI')[1].strip().split()[1].strip()

        for row in self.atoms:
            self.total_charge_penalty_score += float(row.split('!')[1].strip())

        for row in self.bond_parameters:
            if 'PENALTY=' in row:
                self.total_bond_penalty_score += float(row.split('PENALTY=')[1].strip())

        for row in self.angle_parameters:
            if 'PENALTY=' in row:
                self.total_angle_penalty_score += float(row.split('PENALTY=')[1].strip())

        for row in self.dihedral_parameters:
            if 'PENALTY=' in row:
                self.total_dihedral_penalty_score += float(row.split('PENALTY=')[1].strip())

        for row in self.improper_dihedral_parameters:
            if 'PENALTY=' in row:
                self.total_improper_penalty_score += float(row.split('PENALTY=')[1].strip())

    def create_cgenff_dataframes(self):


        '''

        Create the CGenFF Dataframe


        Example Stream output file

        BONDS
        CG2O2  CG312   200.00     1.5220 ! ***** , from CG2O2 CG321, penalty= 55

        ANGLES
        CG312  CG2O2  OG2D1    70.00    125.00   20.00   2.44200 ! ***** , from CG321 CG2O2 OG2D1, penalty= 12

        DIHEDRALS
        OG2D1  CG2O2  CG312  CG312      0.0500  6   180.00 ! ***** , from OG2D1 CG2O2 CG321 CG321, penalty= 67

        '''

        bond_parameter_types = []
        angle_parameter_types = []
        dihedral_parameter_types = []
        improper_parameter_types = []

        bonded_force_constants = []
        bonded_lengths = []

        angle_force_constants = []
        angle_degrees = []

        dihedral_force_constants = []
        dihedral_degrees = []
        dihedral_multipliticies = []

        improper_force_constants = []
        improper_degrees = []
        improper_multipliticies = []

        charge_values = []
        atom_rank_ordering = []
        nonbonded_parameter_types = []

        for row in self.atoms:

            atom_rank = row.split()[1]
            atom_type = row.split()[2]
            charge_value = row.split()[3]

            if len(atom_rank) == 1:

                atom_rank = atom_rank + str(self.counter)

                self.counter += 1

            atom_rank_ordering.append(atom_rank)
            charge_values.append(charge_value)
            nonbonded_parameter_types.append(atom_type)

        for row in self.bond_parameters:

            bond_type = '-'.join(row.split()[0:2])
            bond_length = row.split()[3]
            bond_force_constant = row.split()[2]

            bond_parameter_types.append(bond_type)
            bonded_lengths.append(bond_length)
            bonded_force_constants.append(bond_force_constant)

        for row in self.angle_parameters:

            angle_type = '-'.join(row.split()[0:3])
            angle_parameter_types.append(angle_type)

            angle_degree = row.split()[4]
            angle_force_constant = row.split()[3]

            angle_degrees.append(angle_degree)
            angle_force_constants.append(angle_force_constant)

        for row in self.dihedral_parameters:

            dihedral_type = '-'.join(row.split()[0:4])
            dihedral_parameter_types.append(dihedral_type)

            dihedral_degree = row.split()[6]
            dihedral_multiplicity = row.split()[5]
            dihedral_force_constant = row.split()[4]

            dihedral_degrees.append(dihedral_degree)
            dihedral_multipliticies.append(dihedral_multiplicity)
            dihedral_force_constants.append(dihedral_force_constant)

        for row in self.improper_dihedral_parameters:

            improper_type = '*' + '-'.join(row.split()[0:4])
            improper_parameter_types.append(improper_type)

            improper_degree = row.split()[6]
            improper_multiplicity = row.split()[5]
            improper_force_constant = row.split()[4]

            improper_degrees.append(improper_degree)
            improper_multipliticies.append(improper_multiplicity)
            improper_force_constants.append(improper_force_constant)


        self.bonded_dataframe['atom_types'] = bond_parameter_types + angle_parameter_types + dihedral_parameter_types + improper_parameter_types
        self.bonded_dataframe['force_constants'] = bonded_force_constants + angle_force_constants + dihedral_force_constants + improper_force_constants
        self.bonded_dataframe['harmonic_value'] = bonded_lengths + angle_degrees + dihedral_degrees + improper_degrees
        self.bonded_dataframe['multiplicities'] = (len(bonded_lengths + angle_degrees) * [0]) + dihedral_multipliticies + improper_multipliticies

        self.nonbonded_dataframe['atom_types'] = nonbonded_parameter_types
        self.nonbonded_dataframe['atom_ranks'] = atom_rank_ordering
        self.nonbonded_dataframe['charge'] = charge_values

        self.new_bonded_dataframe = self.bonded_dataframe.copy()
        self.new_nonbonded_dataframe = self.nonbonded_dataframe.copy()

    def update_bonded_dataframe(self, key, new_value, force_constant=True, harmonic_value=False, multiplicity=False):

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

    def update_nonbonded_dataframe(self, key, new_value, charge=True):

        '''

        Update a value in the forcefield dataframes

        Arguments:

            key (String): Bonded key you would like to change
            new_value (String): new value to update
            replace_old_frame (Bool): Whether to replace the old frame
            charge (Bool): Whether the user wants the charge value updated

        '''

        for i, row in self.new_nonbonded_dataframe.iterrows():

            if row['atom_types'] == key:

                if charge:

                    self.new_bonded_dataframe.at[i, 'charge' ] = new_value

    def write_stream_file(self, file_name = 'cgenff_modified.str'):

        '''

        Write the CGenFF Stream File

        Arguments:
            file_name (String): File name of the cgenff stream file

        '''

        out_file = open(file_name, 'w')

        cgenff_header = self._get_cgenff_header()
        cgenff_template = self._get_cgenff_param_template()
        cgenff_end = self._get_cgenff_end()

        # Header

        out_file.write(cgenff_header)

        # Atom Ranks and Charges

        for i, row in self.nonbonded_dataframe.iterrows():

            out_file.write(f'ATOM {row["atom_ranks"]}     {row["atom_types"]}    {row["charge"]}\n' )

        out_file.write('\n')

        # Bond Connectivity

        for row in self.bond_connectivity:
            out_file.write(row + '\n')

        for row in self.improper_atoms:
            out_file.write(row + '\n')

        # CGenFF Template

        out_file.write(cgenff_template)

        out_file.write('\n')
        out_file.write('\n')

        # Bonds

        out_file.write("BONDS")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 2:
                out_file.write(f'{bond_parameter[0]}  {bond_parameter[1]}   {row["force_constants"]}     {row["harmonic_value"]}\n' )

        # Angles

        out_file.write('\n')
        out_file.write("ANGLES")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 3:
                out_file.write(f'{bond_parameter[0]}  {bond_parameter[1]}  {bond_parameter[2]}   {row["harmonic_value"]}     {row["force_constants"]}\n' )


        # DIHEDRALS

        out_file.write('\n')
        out_file.write("DIHEDRALS")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 4 and '*' not in row['atom_types']:
                out_file.write(f'{bond_parameter[0]}  {bond_parameter[1]}  {bond_parameter[2]} {bond_parameter[3]}   {row["force_constants"]}  {row["multiplicities"]}   {row["harmonic_value"]} \n' )

        # IMPROPER DIHEDRALS

        out_file.write('\n')
        out_file.write("IMPROPERS")
        out_file.write('\n')

        for i, row in self.bonded_dataframe.iterrows():

            bond_parameter = row['atom_types'].split('-')

            if len(bond_parameter) == 4 and '*' in row['atom_types']:
                out_file.write(f'{bond_parameter[0].replace("*", "")}  {bond_parameter[1]}  {bond_parameter[2]} {bond_parameter[3]}   {row["force_constants"]}  {row["multiplicities"]}   {row["harmonic_value"]} \n' )

        out_file.write('\n')
        out_file.write(cgenff_end)

    def _get_cgenff_header(self):

        '''

        Return the CGenFF Header

        '''

        cgenff_header = '''
            * Toppar stream file generated by
            * CHARMM General Force Field (CGenFF) program version 2.5
            * For use with CGenFF version 4.5
            *
            
            read rtf card append
            * Topologies generated by
            * CHARMM General Force Field (CGenFF) program version 2.5
            *
            36 1
            
            ! "penalty" is the highest penalty score of the associated parameters.
            ! Penalties lower than 10 indicate the analogy is fair; penalties between 10
            ! and 50 mean some basic validation is recommended; penalties higher than
            ! 50 indicate poor analogy and mandate extensive validation/optimization.
            
            RESI %s          %s
            GROUP
            
            ''' % (self.residue_name, self.molecule_charge)

        return textwrap.dedent(cgenff_header)

    def _get_cgenff_param_template(self):

        '''

        Return the CGenFF Param Card Append.

        '''

        cgenff_template = ''' \
        
        END
        
        read param card flex append
        * Parameters generated by analogy by
        * CHARMM General Force Field (CGenFF) program version 2.5
        *
        
        ! Penalties lower than 10 indicate the analogy is fair; penalties between 10
        ! and 50 mean some basic validation is recommended; penalties higher than
        ! 50 indicate poor analogy and mandate extensive validation/optimization.
        
        '''

        return textwrap.dedent(cgenff_template)

    def _get_cgenff_end(self):


        '''

        Tail end of the CGenFF

        '''

        cgenff_end = '''
        
            END
            
            RETURN
            
        '''

        return textwrap.dedent(cgenff_end).strip()

    def convert_to_rdkit(self):

        '''

        Convert the CGenFF Molecule in RDKit

        '''

        # Blank Molecule

        molecule = Chem.Mol()
        editable_molecule = Chem.EditableMol(molecule)

        # Create an Index Mapping

        atom_index_mapping = {}

        for i, row in enumerate(self.atoms):

            element = row.split()[1]
            atom_index_mapping[element] = i

        for key in list(atom_index_mapping.keys()):

            periodic_table = GetPeriodicTable()
            number = periodic_table.GetAtomicNumber(key[0])
            editable_molecule.AddAtom(Chem.Atom(number))

        for row in self.bond_connectivity:

            atom_1 = row.split()[1]
            atom_2 = row.split()[2]

            editable_molecule.AddBond(
                atom_index_mapping[ atom_1 ],
                atom_index_mapping[ atom_2 ]
            )

        molecule = editable_molecule.GetMol()

        molecule.UpdatePropertyCache()
        Chem.SanitizeMol(molecule)

        AllChem.EmbedMolecule(molecule)
        AllChem.MMFFOptimizeMolecule(molecule)

        return molecule

    def convert_to_smiles(self):

        '''

        Convert the CGenFF to SMILES

        '''

        molecule = self.convert_to_rdkit()

        return Chem.MolToSmiles(molecule)

    def __str__(self):

        return str(self.new_nonbonded_dataframe) + str(self.new_bonded_dataframe)