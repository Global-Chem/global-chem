#!/usr/bin/env python3
#
# GlobalChemExtensions - Drug Design Filters
#
# ------------------------------------------

# RDKit Imports
# -------------

from rdkit import Chem
import rdkit.Chem.Descriptors as Descriptors

class DrugDesignFilters(object):

    __version__ = '0.0.1'

    def __init__(self,
                 smiles_list,
                 lipinski_rule_of_5=False,
                 ghose=False,
                 veber=False,
                 rule_of_3=False,
                 reos=False,
                 drug_like=False,
                 pass_all_filters=False
        ):

        self.smiles_list = smiles_list
        self.molecules = [ Chem.MolFromSmiles(i) for i in self.smiles_list ]

        self.lipinski_rule_of_5 = lipinski_rule_of_5
        self.ghose = ghose
        self.veber = veber
        self.rule_of_3 = rule_of_3
        self.reos = reos
        self.drug_like = drug_like
        self.pass_all_filters = pass_all_filters

    def fetch_attributes(self, molecule):

        '''

        Fetch the Attributes of the molecule

        '''

        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
        number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(molecule)
        molar_refractivity = Chem.Crippen.MolMR(molecule)
        topological_surface_area_mapping = Chem.QED.properties(molecule).PSA
        formal_charge = Chem.rdmolops.GetFormalCharge(molecule)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(molecule)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(molecule)

        return (
            molecular_weight,
            logp,
            h_bond_donor,
            h_bond_acceptors,
            rotatable_bonds,
            number_of_atoms,
            molar_refractivity,
            topological_surface_area_mapping,
            formal_charge,
            heavy_atoms,
            num_of_rings
        )

    def filter(self):

        '''

        Initialize the filtering

        '''

        results = []

        for i in self.molecules:

            (
                molecular_weight,
                logp,
                h_bond_donor,
                h_bond_acceptors,
                rotatable_bonds,
                number_of_atoms,
                molar_refractivity,
                topological_surface_area_mapping,
                formal_charge,
                heavy_atoms,
                num_of_rings
            ) = self.fetch_attributes(i)

            if self.lipinski_rule_of_5 or self.pass_all_filters:

                if molecular_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10and rotatable_bonds <= 5:

                    passed_lipinski = True
                    results.append(Chem.MolToSmiles(i))

            if self.ghose or self.pass_all_filters:

                if molecular_weight >= 160 and molecular_weight <= 480 and logp >= -0.4 and logp <= 5.6 and number_of_atoms >= 20 and number_of_atoms <= 70 and molar_refractivity >= 40 and molar_refractivity <= 130:

                    results.append(Chem.MolToSmiles(i))

            if self.veber or self.pass_all_filters:

                if rotatable_bonds <= 10 and topological_surface_area_mapping <= 140:

                    results.append(Chem.MolToSmiles(i))

            if self.rule_of_3 or self.pass_all_filters:

                if molecular_weight <= 300 and logp <= 3 and h_bond_donor <= 3 and h_bond_acceptors <= 3 and rotatable_bonds <= 3:

                    results.append(Chem.MolToSmiles(i))

            if self.reos or self.pass_all_filters:

                if molecular_weight >= 200 and molecular_weight <= 500 and logp >= int(0 - 5) and logp <= 5 and h_bond_donor >= 0 and h_bond_donor <= 5 and h_bond_acceptors >= 0 and h_bond_acceptors <= 10 and formal_charge >= int(0-2) and formal_charge <= 2 and rotatable_bonds >= 0 and rotatable_bonds <= 8 and heavy_atoms >= 15 and heavy_atoms <= 50:

                    results.append(Chem.MolToSmiles(i))

            if self.drug_like or self.pass_all_filters:

                if molecular_weight < 400 and num_of_rings > 0 and rotatable_bonds < 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and logp < 5:

                    results.append(Chem.MolToSmiles(i))

        return results


