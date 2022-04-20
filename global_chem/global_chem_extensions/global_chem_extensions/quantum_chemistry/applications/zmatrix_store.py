#!/usr/bin/env python3
#
# GlobalChemExtensions - ZMatrix Store

# ------------------------------------

class ZMatrixStore(object):

    '''

    A ZMatrix Store that will house molecules for people to keep downloading if they need it. Will put here in extensions
    to allow for the any conversions later on.

    '''

    __MOLECULES__ = [
        'water',
        'carbon_monoxide',
        'alkyne',
        'cyclopropane',
        'benzene',
    ]

    def __init__(self):

        self.store = {}
        self.build_coordinate_store()

    def get_molecule(self, molecule):

        '''

        Fetch the Molecule

        '''

        if molecule not in self.__MOLECULES__:
            return None
        else:
            return self.store[molecule]

    def build_coordinate_store(self):

        '''

        Builds the Coordinate Store

        '''

        water = """\
        O
        H 1 1.08
        H 1 1.08 2 107.5
        """

        carbon_monoxide = """\
        C
        O 1 1.128
        """

        alkyne = """\
        C          0.99383        0.02433       -0.09699
        C          2.19414        0.02433       -0.09699
        H         -0.07195        0.02433       -0.09699
        H          3.25992        0.02433       -0.09699
        """

        cyclopropane = """\
        C         -0.02818        0.86300        0.05370
        C         -0.73719       -0.45561       -0.06956
        C          0.76431       -0.40718       -0.06957
        H         -0.04424        1.36121        1.01652
        H         -0.04983        1.53422       -0.79697
        H         -1.23857       -0.67874       -1.00439
        H         -1.23417       -0.84996        0.80931
        H          1.27901       -0.59757       -1.00440
        H          1.28568       -0.76871        0.80929
        """

        benzene_xyz = """\
        C         -0.76001        1.16914       -0.00051
        C          0.63286        1.24469       -0.00118
        C          1.39473        0.07653        0.00041
        C          0.76406       -1.16766        0.00273
        C         -0.62881       -1.24322        0.00014
        C         -1.39068       -0.07506       -0.00152
        H         -1.35358        2.07925        0.00053
        H          1.12430        2.21404       -0.00281
        H          2.47994        0.13547       -0.00003
        H          1.35763       -2.07777        0.00630
        H         -1.12024       -2.21257       -0.00047
        H         -2.47589       -0.13400       -0.00346
        """
        benzene = """\
        H
        C 1 a2bd
        C 2 a3bd 1 a3angle
        C 3 a4bd 2 a4angle 1 a4dihedral
        H 4 a5bd 3 a5angle 2 a5dihedral
        C 4 a6bd 3 a6angle 2 a6dihedral
        H 6 a7bd 4 a7angle 3 a7dihedral
        C 6 a8bd 4 a8angle 3 a8dihedral
        H 8 a9bd 6 a9angle 4 a9dihedral
        C 2 a10bd 3 a10angle 4 a10dihedral
        H 10 a11bd 2 a11angle 1 a11dihedral
        H 3 a12bd 2 a12angle 1 a12dihedral
        a2bd = 1.0853
        a3bd = 1.3800
        a4bd = 1.3787
        a5bd = 1.0816
        a6bd = 1.3795
        a7bd = 1.0844
        a8bd = 1.3872
        a9bd = 1.0847
        a10bd = 1.3812
        a11bd = 1.0838
        a12bd = 1.0818
        a3angle = 119.8001
        a4angle = 121.2010
        a5angle = 121.0431
        a6angle = 121.7007
        a7angle = 122.4276
        a8angle = 116.6791
        a9angle = 118.3291
        a10angle = 117.8148
        a11angle = 117.3613
        a12angle = 117.2942
        a4dihedral = -180
        a5dihedral = -180
        a6dihedral = 0
        a7dihedral = 180
        a8dihedral = 0
        a9dihedral = 180
        a10dihedral = 0
        a11dihedral = 0
        a12dihedral = 0
        """

        molecules = [
            water,
            carbon_monoxide,
            alkyne,
            cyclopropane,
            benzene,
        ]

        for i, molecule in enumerate(molecules):

            self.store[self.__MOLECULES__[i]] = molecule

    def get_store(self):

        '''

        Returns the ZMatrix Store

        '''

        return self.store