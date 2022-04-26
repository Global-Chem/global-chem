# GlobalChem Protein

GlobalChem has the concept of a protein into a 1-Dimensional format that can be used later on for protein-smiles to protein-ligand interaction of lexical strings. This way we can perhaps make casual inference between sets of ligands and sets of proteins.&#x20;

**Import the Package**

```
from global_chem_extensions import GlobalChemExtensions

bioinformatics = GlobalChemExtensions().bioinformatics()
```

**Initialize a GlobalChem Protein**

You can initialize a `GlobalChem Protein` with a variety of ways, fetching a pdb  or passing in a `pdb` file as as available in the `Biopandas` package or just passing in a peptide sequence.&#x20;

{% embed url="http://rasbt.github.io/biopandas" %}

{% tabs %}
{% tab title="Code" %}
```
gc_protein = bioinformatics.initialize_globalchem_protein(
        peptide_sequence='AAAA',
)
```
{% endtab %}

{% tab title="Output" %}
```
<PandasPDB> Object
```
{% endtab %}
{% endtabs %}

**Convert to SMILES**

Convert a protein to a smiles string where you can mark the carbonyl carbon or the nitrogen backbone as ionizable centers atom class of `1` or `2.`

{% tabs %}
{% tab title="Code" %}
```
smiles_protein = gc_protein.convert_to_smiles(
    mark_nitrogen_backbone=False,
    mark_carbonyl_carbon_backbone = False,
)
print (smiles_protein)
```
{% endtab %}

{% tab title="Output" %}
```
[N:1]C(C)C([N:1]C(C)C([N:1]C(C)C([N:1]C(C)C([N:1]CC(O)=O)=O)=O)=O)=O
```
{% endtab %}
{% endtabs %}

**Convert to SMARTS**

Convert a protein to a smarts string where you can mark the carbonyl carbon or the nitrogen backbone as atom class of `1` or `2.`

{% tabs %}
{% tab title="Code" %}
```
smarts_protein = gc_protein.convert_to_smarts(
    mark_nitrogen_backbone=False,
    mark_carbonyl_carbon_backbone = False,
)
print (smarts_protein)
```
{% endtab %}

{% tab title="Output" %}
```
[#7:1]-[#6](-[#6])-[#6](-[#7:1]-[#6](-[#6])-[#6](-[#7:1]-[#6](-[#6])-[#6](-[#7:1]-[#6](-[#6])-[#6](-[#7:1]-[#6]-[#6](-[#8])=[#8])=[#8])=[#8])=[#8])=[#8]
```
{% endtab %}
{% endtabs %}

**Convert Ligand to SMILES**

If the PDB has a ligand attached we convert the ligand to SMILES:

{% tabs %}
{% tab title="Code" %}
```
smiles_molecule = gc_protein.convert_ligand_to_smiles()
```
{% endtab %}

{% tab title="Output" %}
```
CS(=O)(=O)N1CCOC(C1)c2csc(n2)c3ccccc3NC(=O)c4[nH]ccn4
```
{% endtab %}
{% endtabs %}

**Filter Criteria to Bostrom's Algorithm**

A useful criteria for medicinal chemists is to separate ligands out of pdb files and determine if the ligand is suitable for the drug binding. The criteria and the code can be found below:

```
1.)  80 < molecular weight (Da) < 750
2.)  10 < number of nonhydrogen atoms < 70
3.)  Must not contain atoms of types other than
     H, C, O, N, F, P, S, Cl, Br, or I
4.)  must contain at least one
     non-carbon/non-hydrogen atom
5.)  must not contain two or more phosphorus atoms
6.)  must not have more than 10 rotatable bonds
7.)  must not be a nucleic acid
8.)  must not be composed only from
     non-lead-like PDB-HET-groupsb
9.)  must not be covalently bound
10.) must not have protein contacts from the crystal
     packing environment in less than 3 Å distance
 11.) must have contacts with protein in less than
      7Å distance
```

{% tabs %}
{% tab title="Code" %}
```
pdb_ids = ['5tc0', '4EHM',]

for pdb in pdb_ids:
    gc_protein = gce.initialize_globalchem_protein(
        pdb_id='5tc0',
    )

    criteria_met = gc_protein.determine_bostrom_ligand_criteria(verbose=True)
    ligand_smiles = gc_protein.ligand_smiles

    print ("Ligand: %s, Criteria Met: %s" % (ligand_smiles, criteria_met))
```
{% endtab %}

{% tab title="Output" %}
```
Passed Check 1 Molecular Weight: 433.08784608800005 
Passed Check 2 Non Hydrogen Atoms: 29 
Passed Check 3 All Atoms Are Within Element Boundaries
Passed Check 4 Non-Hydrogen & Non-Carbon Atoms Present: 11
Passed Check 5 Less than Two Phosphorous atoms Present
Passed Check 6 Less than Ten Rotatable Bonds Present
Passed Check 7 No Nucleic Acid Template Found
Passed Check 9 Not a Covalent Inhibitor
Check 11 Has Contacts within 7 Angstroms
Ligand: CS(=O)(=O)N1CCOC(C1)c2csc(n2)c3ccccc3NC(=O)c4[nH]ccn4, Criteria Met: True
Passed Check 1 Molecular Weight: 358.0688674079999 
Passed Check 2 Non Hydrogen Atoms: 26 
Passed Check 3 All Atoms Are Within Element Boundaries
Passed Check 4 Non-Hydrogen & Non-Carbon Atoms Present: 8
Passed Check 5 Less than Two Phosphorous atoms Present
Passed Check 6 Less than Ten Rotatable Bonds Present
Passed Check 7 No Nucleic Acid Template Found
Failed Check 9
Passed Check 1 Molecular Weight: 358.0688674079999 
Passed Check 2 Non Hydrogen Atoms: 26 
Passed Check 3 All Atoms Are Within Element Boundaries
Passed Check 4 Non-Hydrogen & Non-Carbon Atoms Present: 8
Passed Check 5 Less than Two Phosphorous atoms Present
Passed Check 6 Less than Ten Rotatable Bonds Present
Passed Check 7 No Nucleic Acid Template Found
Failed Check 9
Ligand: Cc1cc(c(c2c1C(=O)Oc3c(c(cc(c3O2)C(=O)O)OC)C)C=O)O, Criteria Met: False
```
{% endtab %}
{% endtabs %}

**Reference**

1.) Boström, Jonas, et al. “Do Structurally Similar Ligands Bind in a Similar Fashion?” Journal of Medicinal Chemistry, vol. 49, no. 23, Nov. 2006, pp. 6716–25.

