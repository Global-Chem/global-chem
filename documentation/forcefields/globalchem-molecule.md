# GlobalChem Molecule



The `GlobalChem` Molecule will be used to house both cheminformatic and forcefield information about a molecule. That way it can serve as a connection piece.&#x20;

**Imports**

```
from global_chem_extensions import GlobalChemExtensions
ff = GlobalChemExtensions().forcefields()
```

**Initialize a GlobalChem Molecule with SMILES**&#x20;

For this example, we are going to use a stream file that is within the global-chem repository as perfluorobutanoic acid.&#x20;

```
global_chem_molecule = ff.initialize_globalchem_molecule(
    'FC(F)(C(F)(C(O)=O)F)C(F)(F)F',
    stream_file='global-chem/example_data/forcefield_parameterization/perfluorobutanoic_acid.str'',
    frcmod_file='gaff2.frcmod',
)
```

{% hint style="info" %}
There are optional parameters to provide the CGenFF stream file or the frcmod file.&#x20;
{% endhint %}

**Determine the Name**

Determine the IUPAC name from SMILES

{% tabs %}
{% tab title="Code" %}
```
global_chem_molecule.determine_name()
name = global_chem_molecule.name
print (name)
```
{% endtab %}

{% tab title="Output" %}
```
'name': 'benzene'
```
{% endtab %}
{% endtabs %}

**Get Attributes**

We take useful cheminformatic descriptors from `RDKit` to retrieve the molecule.&#x20;

{% tabs %}
{% tab title="Code" %}
```
attributes = global_chem_molecule.get_attributes()
for k, v in attributes.items():
  print (f'{k}: {v}')
```
{% endtab %}

{% tab title="Output" %}
```
Attributes: {
'name': 'cyclopentane', 'smiles': 'C1CCCC1',
'molecular_weight': 70.07825032, 
'logp': 1.9505000000000001,
'h_bond_donor': 0,
'h_bond_acceptors': 0,
'rotatable_bonds': 0, 
'number_of_atoms': 5,
'molar_refractivity': 23.084999999999994,
'topological_surface_area_mapping': 0.0,
'formal_charge': 0,
'heavy_atoms': 5,
'num_of_rings': 1
}

```
{% endtab %}
{% endtabs %}

**Convert to Different Interoperable mol or SMILES objects**

To allow for any interoperability with different software we can provide a conversion into their respective objects.&#x20;

{% tabs %}
{% tab title="Code" %}
```
pysmiles_mol = global_chem_molecule.get_pysmiles()
rdkit_mol = global_chem_molecule.get_rdkit_molecule()
partial_smiles_mol = global_chem_molecule.get_partial_smiles()
deep_smiles_mol = global_chem_molecule.encode_deep_smiles()
selfies_mol = global_chem_molecule.encode_selfies()
validate_smiles = global_chem_molecule.validate_molvs()

print ("RDKit Mol: %s" % rdkit_mol)
print ("PySMILES: %s" % pysmiles_mol)
print ("Partial SMILES: %s" % partial_smiles_mol)
print ("DeepSMILES: %s" % deep_smiles_mol)
print ("Selfies Mol: %s" % selfies_mol)
print ("Validation: %s" % validate_smiles)
```
{% endtab %}

{% tab title="Output" %}
```
RDKit Mol: <rdkit.Chem.rdchem.Mol object at 0x7f7259a7cb20>
PySMILES: Graph with 5 nodes and 5 edges
Partial SMILES: Molecule(atoms=['Atom(elem=6,chg=0,idx=0)', 'Atom(elem=6,chg=0,idx=1)', 'Atom(elem=6,chg=0,idx=2)', 'Atom(elem=6,chg=0,idx=3)', 'Atom(elem=6,chg=0,idx=4)'])
DeepSMILES: CCCCC5
Selfies Mol: [C][C][C][C][C][Ring1][Branch1]
Validation: ['INFO: [FragmentValidation] pentane is present']
```
{% endtab %}
{% endtabs %}

{% hint style="info" %}
If a failure at the conversion happens, it must be because perhaps the SMILES is invalid. Please check with the validation module.&#x20;
{% endhint %}

**Draw CGenFF Molecule**

{% tabs %}
{% tab title="Code" %}
```
global_chem_molecule.draw_cgenff_molecule(
    height=1000, 
    width=1000
 )
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/Screen Shot 2022-04-07 at 11.10.13 PM.png>)
{% endtab %}
{% endtabs %}

**Get CGenFF CXSMILES**

CXSMILES is a atoms in the SMILES string followed by a `|` and then a set of features, or metadata, to accommodate each atom. This allows the user to embed whatever information they wish into the molecule.&#x20;

We decided to embed atom-type information since the atom-type ordering is preserved.

{% tabs %}
{% tab title="Code" %}
```
print (global_chem_molecule.get_cgenff_cxsmiles())
```
{% endtab %}

{% tab title="Output" %}
```
O=C(O)C(F)(F)C(F)(F)C(F)(F)F |atomProp:0.atom_type.OG2D1:0.atom_idx.7:1.atom_type.CG2O2:1.atom_idx.5:2.atom_type.OG311:2.atom_idx.6:3.atom_type.CG312:3.atom_idx.3:4.atom_type.FGA2:4.atom_idx.4:5.atom_type.FGA2:5.atom_idx.8:6.atom_type.CG312:6.atom_idx.1:7.atom_type.FGA2:7.atom_idx.0:8.atom_type.FGA2:8.atom_idx.2:9.atom_type.CG302:9.atom_idx.9:10.atom_type.FGA3:10.atom_idx.10:11.atom_type.FGA3:11.atom_idx.11:12.atom_type.FGA3:12.atom_idx.12|
```
{% endtab %}
{% endtabs %}
