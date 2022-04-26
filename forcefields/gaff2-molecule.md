# GAFF2 Molecule

The GAFF2 Molecule follows similar paradigms as the CGenFF Molecule in terms of parsing. GAFF2 also includes a penalty score which makes this more interesting to capture.&#x20;

![](<../.gitbook/assets/Screen Shot 2022-04-03 at 8.18.30 AM.png>)

The equation differs slightly in terms of the harmonics and recent additions have been the impropers. The parser was built to handle _frcmod_ files into dataframes and update their values accordingly.&#x20;

An example of the format can be found here in terms of how the frcmod is defined and what all the numbers mean:

```
BOND
atom type | force constant | bond length
ca-ns  315.20   1.412       same as ca- n, penalty score=  0.0

ANGLE
atom type | force constant | angle degree
c -ns-ca   65.700     123.710   same as c -n -ca, penalty score=  0.0

DIHE
atom type | idivf1 | force constant | phase | periodicity
o -c -ns-ca   4   10.000       180.000           2.000      same as X -c -n -X , penalty score=  0.0

IMPROPER
atom type | force constant | phase | periodicity
ca-ca-ca-ns         1.1          180.0         2.0          Using the default value

NONBON
atom type | rmin half | epsilon
```

**Import the Package**

```
from global_chem_extensions import GlobalChemExtensions
gce = GlobalChemExtensions()
```

**Read in GAFF2 frcmod File**

{% tabs %}
{% tab title="Code" %}
```
gaff2_molecule = gce.initialize_gaff2_molecule('file.frcmod')
print (gaff2_molecule)
```
{% endtab %}
{% endtabs %}

&#x20;**Access Non-Bonded and Bonded DataFrames**

{% tabs %}
{% tab title="Code" %}
```
gaff2_molecule = gce.initialize_gaff2_molecule('file.frcmod')
print (gaff2_molecule.nonbonded_dataframe)
print (gaff2_molecule.bonded_dataframe)
```
{% endtab %}
{% endtabs %}

&#x20;**Update a value in Bonded DataFrames**

Harmonic values refers to anything that is either a bond-length, angle, dihedral or improper angle. Multiplicity only refers to dihedrals or impropers.&#x20;

{% tabs %}
{% tab title="Code" %}
```
gaff2_molecule = gce.initialize_gaff2_molecule('file.frcmod')

gaff2_molecule.update_bonded_dataframe(
    'ca-ca-ns-hn',
     200,
     force_constant=True,
     harmonic_value=False,
     multiplicity=False,
     fluctuation=False,
     )

print (gaff2_molecule.new_bonded_dataframe)
```
{% endtab %}
{% endtabs %}

**Access Manipulated Non-Bonded and Bonded DataFrames**

The values to update are the Rmin and the epsilon values.

{% tabs %}
{% tab title="Code" %}
```
gaff2_molecule = gce.initialize_gaff2_molecule('file.frcmod')

gaff2_molecule.update_nonbonded_dataframe(
'ca-ca',
 3,
 rmin=True,
 epsilon=False,
 )
print (gaff2_molecule.new_nonbonded_dataframe)

```
{% endtab %}
{% endtabs %}

**Write a new Frcmod GAFF2 File**

{% tabs %}
{% tab title="Code" %}
```
gaff2_molecule = gce.initialize_gaff2_molecule('file.frcmod')
cgenff_molecule.write_frcmod_file()
```
{% endtab %}

{% tab title="Output" %}
```
Remark line goes here
MASS
BONDS
ca-ns	1.412	315.20
c-ns	1.379	356.20
hn-ns	1.013	527.30

ANGLE
c-ns-ca	123.710	65.700
ca-ns-hn	116.000	48.000
ns-c-o	123.050	113.800
c-ns-hn	117.550	48.700
ca-ca-ns	120.190	85.600
c3-c-ns	115.180	84.300

DIHE
o-c-ns-ca	4	10.000	180.000	2.000
c3-c-ns-ca	1	0.750	180.000	-2.000
c3-c-ns-ca	1	0.500	0.000	3.000
o-c-ns-hn	1	2.500	180.000	-2.000
o-c-ns-hn	1	2.000	0.000	1.000
ca-ca-ns-c	1	0.950	180.000	2.000
ns-c-c3-c3	1	0.000	180.000	-4.000
ns-c-c3-c3	1	0.710	180.000	2.000
ca-ca-ns-hn	4	1.800	180.000	2.000
c3-c-ns-hn	4	10.000	180.000	2.000

IMPROPER
ca-ca-ca-ns-1.1	180.0	2.0 
ca-ca-ca-ha-1.1	180.0	2.0 
c3-ns-c-o-10.5	180.0	2.0 
c-ca-ns-hn-1.1	180.0	2.0 
ca-h4-ca-nb-1.1	180.0	2.0 


NONBON
ca-h4	1.1	180.0
ca-h4	1.1	180.0
```
{% endtab %}
{% endtabs %}

**Accessing Properties on GAFF2 Molecule**

```
print (gaff2_molecule.bond_parameters)
print (gaff2_molecule.angle_parameters)
print (gaff2_molecule.dihedral_parameters)
print (gaff2_molecule.nonbonded_parameters)
print (gaff2_molecule.improper_parameters)
```

