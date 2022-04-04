#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Molecule

# -------------------------------------------

# Imports
# -------

from rdkit import Chem
import rdkit.Chem.Descriptors as Descriptors

class GlobalChemMolecule(object):


    __version__ = '0.0.1'

    def __init__(self, smiles):

        '''

        Arguments:
            smiles (String): Chem.MolFromSmiles(i)

        '''

        self.smiles = smiles
        self.molecule = Chem.MolFromSmiles(self.smiles)

        self.attributes = {}
        self._determine_attributes()

    def get_attributes(self):

        '''

        Returns the attributes

        '''

        return self.attributes

    def _determine_attributes(self):

        '''

        Fetch the Attributes of the molecule

        '''


        molecular_weight = Descriptors.ExactMolWt(self.molecule)
        logp = Descriptors.MolLogP(self.molecule)
        h_bond_donor = Descriptors.NumHDonors(self.molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(self.molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(self.molecule)
        number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(self.molecule)
        molar_refractivity = Chem.Crippen.MolMR(self.molecule)
        topological_surface_area_mapping = Chem.QED.properties(self.molecule).PSA
        formal_charge = Chem.rdmolops.GetFormalCharge(self.molecule)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(self.molecule)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(self.molecule)


        self.attributes['molecular_weight'] = molecular_weight
        self.attributes['logp'] = logp
        self.attributes['h_bond_donor'] = h_bond_donor
        self.attributes['h_bond_acceptors'] = h_bond_acceptors
        self.attributes['rotatable_bonds'] = rotatable_bonds
        self.attributes['number_of_atoms'] = number_of_atoms
        self.attributes['molar_refractivity'] = molar_refractivity
        self.attributes['topological_surface_area_mapping'] = topological_surface_area_mapping
        self.attributes['formal_charge'] = formal_charge
        self.attributes['heavy_atoms'] = heavy_atoms
        self.attributes['num_of_rings'] = num_of_rings

