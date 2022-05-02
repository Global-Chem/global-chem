#!/usr/bin/env python3
#
# GlobalChemExtensions - Psi4 Adapter
#
# -----------------------------------

# Imports
# -------

import os
import moly

# Quantum Chemistry Imports
# -------------------------

import psi4

class Psi4Parser(object):

    '''

    Parser Object for the Psi4 Output File

    '''

    def __init__(self,
                 input_file = None,
                 cube_files = None,
                 wave_function = None,
                 ):

        # Attributes

        if input_file is not None:
            self.input_file = open(input_file, 'r').read()
            self.lines = self.input_file.split('\n')

        self.cube_files = cube_files

        # Energy Terms

        self.nuclear_repulsion_energy = None
        self.one_electron_energy = None
        self.two_electron_energy = None
        self.total_energy = None

        # Dipole Terms

        self.nuclear_dipole_moment = None
        self.electronic_dipole_moment = None
        self.dipole_moment_atomic_units = None
        self.dipole_moment_debye = None

        # Frequency Terms

        self.frequency_all_modes = None

        # Wave Function

        self.wave_function = wave_function

    def get_energy_contributions(self):

        '''

        Variables Set:

            self.nuclear_dipole_moment,
            self.electronic_dipole_moment,
            self.dipole_moment_atomic_units,
            self.dipole_moment_debye

        Returns:

            dipoles (Dict): Dipole moments and their respective x, y, z vectors.

        Returns the energy contributions for each term

        '''

        # Will iterate and get the last value

        for line in self.lines:

            if 'Nuclear Repulsion Energy =' in line:
                self.nuclear_repulsion_energy = float(line.strip().split()[-1])

            if 'One-Electron Energy =' in line:
                self.one_electron_energy = float(line.strip().split()[-1])

            if 'Two-Electron Energy =' in line:
                self.two_electron_energy = float(line.strip().split()[-1])

            if 'Total Energy =' in line:
                self.total_energy = float(line.strip().split()[-1])

        energies = {
            'nuclear_repulsion_energy': self.nuclear_repulsion_energy,
            'one_electron_energy': self.one_electron_energy,
            'two_electron_energy': self.two_electron_energy,
            'total_energy': self.total_energy
        }

        return energies

    def get_dipole_contributions(self):

        '''

        Get the Dipole moments

        Variables Set:

            self.nuclear_dipole_moment,
            self.electronic_dipole_moment,
            self.dipole_moment_atomic_units,
            self.dipole_moment_debye

        Returns:

            dipoles (Dict): Dipole moments and their respective x, y, z vectors.

        '''

        # Will iterate and get the last value

        for i, line in enumerate(self.lines):


            if 'Nuclear Dipole Moment: [e a0]' in line:

                line = self.lines[i + 1].strip()

                self.nuclear_dipole_moment = [
                    line.split('X:')[1].split()[0],
                    line.split('Y:')[1].split()[0],
                    line.split('Z:')[1].split()[0]
                ]

            if 'Electronic Dipole Moment: [e a0]' in line:

                line = self.lines[i + 1].strip()

                self.electronic_dipole_moment = [
                    line.split('X:')[1].split()[0],
                    line.split('Y:')[1].split()[0],
                    line.split('Z:')[1].split()[0]
                ]

            if 'Dipole Moment: [e a0]' in line:

                line = self.lines[i + 1].strip()

                self.dipole_moment_atomic_units = [
                    line.split('X:')[1].split()[0],
                    line.split('Y:')[1].split()[0],
                    line.split('Z:')[1].split()[0],
                    line.split('Total:')[1].split()[0],
                ]

            if 'Dipole Moment: [D]' in line:

                line = self.lines[i + 1].strip()

                self.dipole_moment_debye = [
                    line.split('X:')[1].split()[0],
                    line.split('Y:')[1].split()[0],
                    line.split('Z:')[1].split()[0],
                    line.split('Total:')[1].split()[0],
                ]


        dipoles = {
            'nuclear_dipole_moment': self.nuclear_dipole_moment,
            'electronic_dipole_moment': self.electronic_dipole_moment,
            'dipole_moment_atomic_units': self.dipole_moment_atomic_units,
            'dipole_moment_debye': self.dipole_moment_debye
        }

        return dipoles

    def get_frequencies(self):

        '''

        Get the Frequencies for each term

        '''

        modes = self.input_file.split('post-proj  all modes:[')[1].split(']')[0].replace("'", "").split()
        self.frequency_all_modes = modes

        return modes

    def moly_plot_molecular_orbital(self, cube_file, name='orbital', opacity=0.1, show=False):

        '''

        Plot the Molecular Orbital

        Arguments:
            cube_file (String): File Path to the Cube File
            name (String): Name of the output file
            show (Bool): whether to show it in the jupyter notebook or save the file
            opacity (Float): opacity of the orbitals

        '''

        fig = moly.Figure()
        fig.add_cube(cube_file,
                     iso=0.03,
                     colorscale="rdbu",
                     opacity=opacity)

        if show:
            fig.show()
        else:
            fig.fig.write_image(os.path.join(
                '%s.png' % name)
            )

    def psi4_cubeprop(
            self,
            path,
            orbitals = [],
            number_occupied = 0,
            number_virtual = 0,
            density=False,
    ):

        """
        Run a psi4 cubeprop computation to generate cube files from a given Wavefunction object
        By default this function plots from the HOMO -2 to the LUMO + 2

        Arguments:

            path (String): file path to the output directory of the cube files
            orbitals (List): Any specific orbitals you want included,
            number_occupied (Int): Number of Molecule Orbitals Occupied by the Electron
            number_virtual (Int): Number of Virtual Orbitals Occupied by the Electron
            density (Bool): Density of the cubeprop
        """

        import os.path

        cubeprop_tasks = []

        if isinstance(orbitals, str):
            if orbitals == 'frontier_orbitals':
                cubeprop_tasks.append('FRONTIER_ORBITALS')
        else:
            cubeprop_tasks.append('ORBITALS')

            if number_occupied + number_virtual > 0:
                alpha_electrons = self.wave_function.nalpha()
                molecule = self.wave_function.nmo()
                min_orb = max(1, alpha_electrons + 1 - number_occupied)
                max_orb = min(molecule, alpha_electrons + number_virtual)
                orbitals = [k for k in range(min_orb, max_orb + 1)]

            print(f'Preparing cube files for orbitals: {", ".join([str(orbital) for orbital in orbitals])}')

        if density:
            cubeprop_tasks.append('DENSITY')

        if not os.path.exists(path):
            os.makedirs(path)

        psi4.set_options(
            {
                'CUBEPROP_TASKS': cubeprop_tasks,
                'CUBEPROP_ORBITALS': orbitals,
                'CUBEPROP_FILEPATH': path
            }
        )

        psi4.cubeprop(self.wave_function)
