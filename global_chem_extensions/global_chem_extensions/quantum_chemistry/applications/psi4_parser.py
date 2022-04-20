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
                 input_file,
                 cube_files = None,
                 ):

        # Attributes

        self.input_file = open(input_file, 'r').read()
        self.cube_files = cube_files

        self.lines = self.input_file.split('\n')

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

            if 'Nuclear Repulsion Energy' in line:
                self.nuclear_repulsion_energy = float(line.strip().split()[-1])

            if 'One-Electron Energy' in line:
                self.one_electron_energy = float(line.strip().split()[-1])

            if 'Two-Electron Energy' in line:
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

    def moly_plot_molecular_orbital(self, out_directory=None, show=False):

        '''

        Plot the Molecular Orbital

        '''

        for orbital in self.cube_files:
            fig = moly.Figure()
            fig.add_cube(orbital, iso=0.03, colorscale="rdbu", opacity=0.2)

            fig.fig.write_image(os.path.join(
                out_directory, '%s.png' % orbital)
            )
            
            if show:
                fig.show()