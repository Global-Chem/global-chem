#!/usr/bin/env python3
#
# Intialization of the package
#
# ----------------------------


# ForceFields

from global_chem_extensions.global_chem_extensions.forcefields.cgenff.cgenff_parser import CGenFFMolecule

# Base Package

from global_chem_extensions.global_chem_extensions import GlobalChemExtensions

__all__ = ['GlobalChemExtensions', 'CGenFFMolecule']

__name__ = 'GlobalChemExtensions'