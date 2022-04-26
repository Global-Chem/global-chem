# CGenFF Molecule



We can parse the CGenFF stream file into a Molecule object within `GlobalChem` this allows for easy parameter manipulation and a possible conversion route into other tools. CGenFF Molecules are based on the CHARMM Potential Energy Function with their respective classifiers stemmed from the equation:

![](<../.gitbook/assets/image (3).png>)

The CGenFF Parser was used to build and manipulate CGenFF output stream files. And also provide a possible route into other connections.&#x20;

Parsing the data is self-explanatory but the conversion into SMILES is still underway. Currently, the algorithm is missing bond-orders since the CGenFF file doesn't retain that information in the output.&#x20;

**Read in CGenFF Stream File**

{% tabs %}
{% tab title="Code" %}
```
from global_chem_extensions import GlobalChemExtensions
ff = GlobalChemExtensions().forcefields()
```
{% endtab %}

{% tab title="Output" %}
```
   atom_types atom_ranks  charge
0        FGA2         F1  -0.190
1       CG312         C1   0.439
2        FGA2         F2  -0.190
3       CG312         C2   0.424
4        FGA2         F3  -0.216
5       CG2O2         C3   0.758
6       OG311         O1  -0.597
7       OG2D1         O2  -0.557
8        FGA2         F4  -0.216
9       CG302         C4   0.335
10       FGA3         F5  -0.140
11       FGA3         F6  -0.140
12       FGA3         F7  -0.140
13       HGP1         H1   0.430                  atom_types force_constants harmonic_value multiplicities
0                CG2O2-CG312          200.00         1.5220              0
1                CG302-CG312          250.00         1.5200              0
2                CG312-CG312          198.00         1.5200              0
3          CG312-CG2O2-OG2D1           70.00         125.00              0
4          CG312-CG2O2-OG311           55.00         110.50              0
5           CG312-CG302-FGA3           42.00         112.00              0
6          CG2O2-CG312-CG312           52.00         108.00              0
7           CG2O2-CG312-FGA2           50.00         115.00              0
8          CG302-CG312-CG312           50.00         112.00              0
9           CG302-CG312-FGA2           50.00         112.00              0
10          CG312-CG312-FGA2           50.00         112.00              0
11   OG2D1-CG2O2-CG312-CG312          0.0500         180.00              6
12    OG2D1-CG2O2-CG312-FGA2          0.1000           0.00              1
13    OG2D1-CG2O2-CG312-FGA2          0.9800         180.00              2
14   OG311-CG2O2-CG312-CG312          0.0000         180.00              6
15    OG311-CG2O2-CG312-FGA2          0.1700         180.00              1
16    OG311-CG2O2-CG312-FGA2          0.4600         180.00              2
17    CG312-CG2O2-OG311-HGP1          2.0500         180.00              2
18    FGA3-CG302-CG312-CG312          0.2500           0.00              3
19     FGA3-CG302-CG312-FGA2          0.2500           0.00              3
20   CG2O2-CG312-CG312-CG302          0.2100         180.00              1
21   CG2O2-CG312-CG312-CG302          0.3900           0.00              2
22   CG2O2-CG312-CG312-CG302          0.3500         180.00              3
23   CG2O2-CG312-CG312-CG302          0.1100           0.00              4
24   CG2O2-CG312-CG312-CG302          0.0900         180.00              6
25    CG2O2-CG312-CG312-FGA2          2.0000         180.00              1
26    CG2O2-CG312-CG312-FGA2          0.8000           0.00              2
27    CG302-CG312-CG312-FGA2          0.2500           0.00              3
28     FGA2-CG312-CG312-FGA2          0.2500           0.00              3
29  *CG2O2-CG312-OG2D1-OG311         65.0000           0.00              0
```
{% endtab %}
{% endtabs %}

&#x20;**Access Non-Bonded and Bonded DataFrames**

{% tabs %}
{% tab title="Code" %}
```
cgenff_molecule = ff.initialize_cgenff_molecule('file.str')
print (cgenff_molecule.nonbonded_dataframe)
print (cgenff_molecule.bonded_dataframe)
```
{% endtab %}

{% tab title="Output" %}
```
   atom_types atom_ranks  charge
0        FGA2         F1  -0.190
1       CG312         C1   0.439
2        FGA2         F2  -0.190
3       CG312         C2   0.424
4        FGA2         F3  -0.216
5       CG2O2         C3   0.758
6       OG311         O1  -0.597
7       OG2D1         O2  -0.557
8        FGA2         F4  -0.216
9       CG302         C4   0.335
10       FGA3         F5  -0.140
11       FGA3         F6  -0.140
12       FGA3         F7  -0.140
13       HGP1         H1   0.430
                  atom_types force_constants harmonic_value multiplicities
0                CG2O2-CG312          200.00         1.5220              0
1                CG302-CG312          250.00         1.5200              0
2                CG312-CG312          198.00         1.5200              0
3          CG312-CG2O2-OG2D1           70.00         125.00              0
4          CG312-CG2O2-OG311           55.00         110.50              0
5           CG312-CG302-FGA3           42.00         112.00              0
6          CG2O2-CG312-CG312           52.00         108.00              0
7           CG2O2-CG312-FGA2           50.00         115.00              0
8          CG302-CG312-CG312           50.00         112.00              0
9           CG302-CG312-FGA2           50.00         112.00              0
10          CG312-CG312-FGA2           50.00         112.00              0
11   OG2D1-CG2O2-CG312-CG312          0.0500         180.00              6
12    OG2D1-CG2O2-CG312-FGA2          0.1000           0.00              1
13    OG2D1-CG2O2-CG312-FGA2          0.9800         180.00              2
14   OG311-CG2O2-CG312-CG312          0.0000         180.00              6
15    OG311-CG2O2-CG312-FGA2          0.1700         180.00              1
16    OG311-CG2O2-CG312-FGA2          0.4600         180.00              2
17    CG312-CG2O2-OG311-HGP1          2.0500         180.00              2
18    FGA3-CG302-CG312-CG312          0.2500           0.00              3
19     FGA3-CG302-CG312-FGA2          0.2500           0.00              3
20   CG2O2-CG312-CG312-CG302          0.2100         180.00              1
21   CG2O2-CG312-CG312-CG302          0.3900           0.00              2
22   CG2O2-CG312-CG312-CG302          0.3500         180.00              3
23   CG2O2-CG312-CG312-CG302          0.1100           0.00              4
24   CG2O2-CG312-CG312-CG302          0.0900         180.00              6
25    CG2O2-CG312-CG312-FGA2          2.0000         180.00              1
26    CG2O2-CG312-CG312-FGA2          0.8000           0.00              2
27    CG302-CG312-CG312-FGA2          0.2500           0.00              3
28     FGA2-CG312-CG312-FGA2          0.2500           0.00              3
29  *CG2O2-CG312-OG2D1-OG311         65.0000           0.00              0

```
{% endtab %}
{% endtabs %}

&#x20;**Update a value in Bonded DataFrames**

Harmonic values refers to anything that is either a bond-length, angle, dihedral or improper angle. Multiplicity only refers to dihedrals or impropers.&#x20;

{% tabs %}
{% tab title="Code" %}
```
cgenff_molecule.update_bonded_dataframe(
    'CG302-CG312-CG312',
     200,
     force_constant=True,
     harmonic_value=False,
     multiplicity=False
     )
print (cgenff_molecule.new_bonded_dataframe)
```
{% endtab %}

{% tab title="Output" %}
```
                  atom_types force_constants harmonic_value multiplicities
0                CG2O2-CG312          200.00         1.5220              0
1                CG302-CG312          250.00         1.5200              0
2                CG312-CG312          198.00         1.5200              0
3          CG312-CG2O2-OG2D1           70.00         125.00              0
4          CG312-CG2O2-OG311           55.00         110.50              0
5           CG312-CG302-FGA3           42.00         112.00              0
6          CG2O2-CG312-CG312           52.00         108.00              0
7           CG2O2-CG312-FGA2           50.00         115.00              0
8          CG302-CG312-CG312             200         112.00              0
9           CG302-CG312-FGA2           50.00         112.00              0
10          CG312-CG312-FGA2           50.00         112.00              0
11   OG2D1-CG2O2-CG312-CG312          0.0500         180.00              6
12    OG2D1-CG2O2-CG312-FGA2          0.1000           0.00              1
13    OG2D1-CG2O2-CG312-FGA2          0.9800         180.00              2
14   OG311-CG2O2-CG312-CG312          0.0000         180.00              6
15    OG311-CG2O2-CG312-FGA2          0.1700         180.00              1
16    OG311-CG2O2-CG312-FGA2          0.4600         180.00              2
17    CG312-CG2O2-OG311-HGP1          2.0500         180.00              2
18    FGA3-CG302-CG312-CG312          0.2500           0.00              3
19     FGA3-CG302-CG312-FGA2          0.2500           0.00              3
20   CG2O2-CG312-CG312-CG302          0.2100         180.00              1
21   CG2O2-CG312-CG312-CG302          0.3900           0.00              2
22   CG2O2-CG312-CG312-CG302          0.3500         180.00              3
23   CG2O2-CG312-CG312-CG302          0.1100           0.00              4
24   CG2O2-CG312-CG312-CG302          0.0900         180.00              6
25    CG2O2-CG312-CG312-FGA2          2.0000         180.00              1
26    CG2O2-CG312-CG312-FGA2          0.8000           0.00              2
27    CG302-CG312-CG312-FGA2          0.2500           0.00              3
28     FGA2-CG312-CG312-FGA2          0.2500           0.00              3
29  *CG2O2-CG312-OG2D1-OG311         65.0000           0.00              0
```
{% endtab %}
{% endtabs %}

**Access Manipulated Non-Bonded and Bonded DataFrames**

The only value to update, for now, is the charge. When expanding into Drude we can adjust the parameter here

{% tabs %}
{% tab title="Code" %}
```
cgenff_molecule.update_nonbonded_dataframe('CG302', -0.2, charge=True)
print (cgenff_molecule.new_nonbonded_dataframe)
```
{% endtab %}

{% tab title="Output" %}
```
   atom_types atom_ranks  charge
0        FGA2         F1  -0.190
1       CG312         C1   0.439
2        FGA2         F2  -0.190
3       CG312         C2   0.424
4        FGA2         F3  -0.216
5       CG2O2         C3   0.758
6       OG311         O1  -0.597
7       OG2D1         O2  -0.557
8        FGA2         F4  -0.216
9       CG302         C4   0.335
10       FGA3         F5  -0.140
11       FGA3         F6  -0.140
12       FGA3         F7  -0.140
13       HGP1         H1   0.430
```
{% endtab %}
{% endtabs %}

**Write a new stream CGenFF File**

{% tabs %}
{% tab title="Code" %}
```
cgenff_molecule.write_stream_file()
```
{% endtab %}

{% tab title="Output" %}
```

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

RESI *****          0.000
GROUP

ATOM F1     FGA2    -0.190
ATOM C1     CG312    0.439
ATOM F2     FGA2    -0.190
ATOM C2     CG312    0.424
ATOM F3     FGA2    -0.216
ATOM C3     CG2O2    0.758
ATOM O1     OG311    -0.597
ATOM O2     OG2D1    -0.557
ATOM F4     FGA2    -0.216
ATOM C4     CG302    0.335
ATOM F5     FGA3    -0.140
ATOM F6     FGA3    -0.140
ATOM F7     FGA3    -0.140
ATOM H1     HGP1    0.430

BOND F1   C1
BOND C1   F2
BOND C1   C2
BOND C2   F3
BOND C2   C3
BOND C3   O1
BOND C3   O2
BOND C2   F4
BOND C1   C4
BOND C4   F5
BOND C4   F6
BOND C4   F7
BOND O1   H
IMPR C3     C2     O2     O1

END

read param card flex append
* Parameters generated by analogy by
* CHARMM General Force Field (CGenFF) program version 2.5
*

! Penalties lower than 10 indicate the analogy is fair; penalties between 10
! and 50 mean some basic validation is recommended; penalties higher than
! 50 indicate poor analogy and mandate extensive validation/optimization.



BONDS
CG2O2  CG312   200.00     1.5220
CG302  CG312   250.00     1.5200
CG312  CG312   198.00     1.5200

ANGLES
CG312  CG2O2  OG2D1   125.00     70.00
CG312  CG2O2  OG311   110.50     55.00
CG312  CG302  FGA3   112.00     42.00
CG2O2  CG312  CG312   108.00     52.00
CG2O2  CG312  FGA2   115.00     50.00
CG302  CG312  CG312   112.00     50.00
CG302  CG312  FGA2   112.00     50.00
CG312  CG312  FGA2   112.00     50.00

DIHEDRALS
OG2D1  CG2O2  CG312 CG312   0.0500  6   180.00 
OG2D1  CG2O2  CG312 FGA2   0.1000  1   0.00 
OG2D1  CG2O2  CG312 FGA2   0.9800  2   180.00 
OG311  CG2O2  CG312 CG312   0.0000  6   180.00 
OG311  CG2O2  CG312 FGA2   0.1700  1   180.00 
OG311  CG2O2  CG312 FGA2   0.4600  2   180.00 
CG312  CG2O2  OG311 HGP1   2.0500  2   180.00 
FGA3  CG302  CG312 CG312   0.2500  3   0.00 
FGA3  CG302  CG312 FGA2   0.2500  3   0.00 
CG2O2  CG312  CG312 CG302   0.2100  1   180.00 
CG2O2  CG312  CG312 CG302   0.3900  2   0.00 
CG2O2  CG312  CG312 CG302   0.3500  3   180.00 
CG2O2  CG312  CG312 CG302   0.1100  4   0.00 
CG2O2  CG312  CG312 CG302   0.0900  6   180.00 
CG2O2  CG312  CG312 FGA2   2.0000  1   180.00 
CG2O2  CG312  CG312 FGA2   0.8000  2   0.00 
CG302  CG312  CG312 FGA2   0.2500  3   0.00 
FGA2  CG312  CG312 FGA2   0.2500  3   0.00 

IMPROPERS
CG2O2  CG312  OG2D1 OG311   65.0000  0   0.00 

END

RETURN
```
{% endtab %}
{% endtabs %}

**Convert to RDKit & SMILES**&#x20;

{% hint style="info" %}
We are missing bond orders if you are using this feature
{% endhint %}

{% tabs %}
{% tab title="Code" %}
```
molecule = cgenff_molecule.convert_to_rdkit()
molecule = cgenff_molecule.convert_to_smiles()
```
{% endtab %}

{% tab title="Output" %}
```
[HH]~O~C(~O)~C(~F)(~F)~C(~F)(~F)~C(~F)(~F)~F
```
{% endtab %}
{% endtabs %}

The conversion is still tricky here since we don't have bond orders. Information needs to come from somewhere and can be hard to guess by atomic valence with no connections.&#x20;

```
molecule = Chem.Mol()
editable_molecule = Chem.EditableMol(molecule)

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
```

**Accessing Properties on CGenFF Molecule**

```
print (cgenff_molecule.bond_parameters)
print (cgenff_molecule.angle_parameters)
print (cgenff_molecule.dihedral_parameters)
print (cgenff_molecule.atom_rows)
print (cgenff_molecule.improper_parameters)
```
