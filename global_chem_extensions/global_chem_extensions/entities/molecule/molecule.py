#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Molecule

# -------------------------------------------

# Imports
# -------

from IPython.display import SVG

# RDKit Imports
# -------------

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import rdkit.Chem.Descriptors as Descriptors

# Local Imports
# -------------

from global_chem import GlobalChem

class GlobalChemMolecule(object):


    __version__ = '0.0.1'

    def __init__(
            self,
            smiles,
            name = None,
            cgenff_molecule = None,
            gaff2_molecule = None
    ):

        '''

        Arguments:
            smiles (String): Chem.MolFromSmiles(i)
            cgenff_stream_file (String): Whether the user would like it with a cgenff entity attached
            gaff_frcmod_file (String): Whether the user would like it with a gaff2 entity attached
            cgenff_molecule (CGenFFMolecule Object): Quick hack to pass it and not deal with file path issues
            gaff2_molecule (GaFF2Molecule Object): Quick hack to pass the Gaff2 Molecule and not deal with file path issues

        '''

        # Class Attributes

        self.name = name
        self.smiles = smiles

        # Cheminformatic Molecules

        self.molecule = Chem.MolFromSmiles(self.smiles)

        # ForceField Molecules

        self.cgenff_molecule = cgenff_molecule
        self.gaff2_molecule = gaff2_molecule

        # Molecule Attributes

        self.attributes = {}

        self.attributes['name'] = self.name
        self.attributes['smiles'] = self.smiles

        self._determine_attributes()

        # Converters

        self.converter = None

    def get_attributes(self):

        '''

        Returns the attributes

        '''

        return self.attributes

    def get_rdkit_molecule(self):

        '''

        Returns the RDKit Molecule

        '''

        return Chem.MolFromSmiles(self.smiles)

    def get_cgenff_molecule(self):

        '''

        Returns:
            cgenff_molecule (GlobalChem Object): the molecule object
        '''

        return self.cgenff_molecule

    def get_gaff2_molecule(self):

        '''

        Returns:
            gaff2_molecule (GlobalChem Object): the molecule object
        '''

        return self.gaff2_molecule

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

    def determine_name(self):

        '''

        Determine the Name of a molecule if available

        '''


        gc = GlobalChem()
        smiles_list = gc.get_all_smiles()
        name_list = gc.get_all_names()

        name_found = False

        for i, smiles in enumerate(smiles_list):

            name = name_list[i]

            if smiles == self.smiles and not name_found:
                self.name = name
                self.attributes['name'] = name
                name_found = True

    def set_name(self, name):

        '''

        Arguments:
            name (String): set the name of the compound

        '''

        self.name = name

    def get_partial_smiles(self, partial=None):

        '''

        Get's the partial SMILES mol object

        Arguments:
            partial (Bool): whether the user is passing in partial SMILES
        '''

        import partialsmiles as ps

        return ps.ParseSmiles(self.smiles, partial=partial)

    def get_pysmiles(self):

        '''

        Gets the PySMILES mol object

        '''

        from pysmiles import read_smiles

        return read_smiles(self.smiles)

    def encode_deep_smiles(self):

        '''

        Gets the DeepSMILES mol Object

        '''

        import deepsmiles

        self.converter = deepsmiles.Converter(rings=True, branches=True)

        return self.converter.encode(self.smiles)

    def decode_deep_smiles(self):

        '''

        Decodes the DeepSMILES molObject

        '''

        import deepsmiles

        self.converter = deepsmiles.Converter(rings=True, branches=True)

        return self.converter.decode(self.smiles)

    def encode_selfies(self):

        '''

        Encodes into the Selfies Object

        '''

        import selfies as sf

        return sf.encoder(self.smiles)

    def decode_selfies(self):

        '''

        Decodes the Selfies Object

        '''

        import selfies as sf

        return sf.decoder(self.smiles)

    def validate_molvs(self):

        '''

        Validate the SMILES with MOLVS

        '''

        from molvs import validate_smiles

        return validate_smiles(self.smiles)

    def draw_cgenff_molecule(self, height=400, width=400):

        '''

        Draw the CGenFF Molecule

        Returns:
            SVG (IPython Object): SVG Object for the jupyter notebook
            molSize (Tuple): Size of the molecule

        '''

        if self.cgenff_molecule:

            hetereo_atom_types = [i.split()[2] for i in self.cgenff_molecule.atoms if i.split()[2][0] != 'H']

            return SVG(GlobalChemMolecule.mol_to_svg(self.molecule, hetereo_atom_types, self.name, molSize=(height, width)))

    @staticmethod
    def mol_to_svg(mol, atom_types, name, molSize = (400,400), kekulize = True):

        '''

        Mol to SVG Function

        '''

        mc = Chem.Mol(mol.ToBinary())

        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)

        for i in range(0, len(atom_types)):

            atom_type = atom_types[i]
            mc.GetAtomWithIdx(i).SetProp("atomNote", atom_type)

        drawer = rdMolDraw2D.MolDraw2DSVG( molSize[0], molSize[1] )
        drawer.drawOptions().annotationFontScale = 0.5
        drawer.drawOptions().centreMoleculesBeforeDrawing = False
        drawer.drawOptions().useBWAtomPalette()

        drawer.DrawMoleculeWithHighlights(mc, name , {}, {}, {}, {})
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        file = open('%s.svg' % name, 'w')
        file.write(svg)
        file.close()

        return svg.replace('svg:','')

    def get_cgenff_cxsmiles(self):

        '''

        Get the CGenFF CX SMILES

        '''

        if self.cgenff_molecule:

            hetereo_atom_types = [i.split()[2] for i in self.cgenff_molecule.atoms if i.split()[2][0] != 'H']

            for idx, atom in enumerate(self.molecule.GetAtoms()):
                atom.SetProp('atom_type', hetereo_atom_types[idx])
                atom.SetIntProp('atom_idx', idx)

        return Chem.MolToCXSmiles(self.molecule)
