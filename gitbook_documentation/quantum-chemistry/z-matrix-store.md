# Z-Matrix Store



It would be handy to have a set of z-matrix molecules that are common for people to use and manipulate. Fair warning this will be removing cartesian coordinates within the dataset really soon.&#x20;

First initialize the Package

```
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChemExtensions().quantum_chemistry()
```

**Fetch The ZMatrix Store Object**

{% tabs %}
{% tab title="Code" %}
```
store_object = gc.get_zmatrix_store()
store_contents = store_object.get_store()
print (store_contents)
```
{% endtab %}
{% endtabs %}

**Get Molecule**

Fetch a molecule's zmatrix.

{% tabs %}
{% tab title="Code" %}
```
store_object.get_molecule('benzene')
```
{% endtab %}

{% tab title="Output" %}
```
'        H\n        C 1 a2bd\n        C 2 a3bd 1 a3angle\n        C 3 a4bd 2 a4angle 1 a4dihedral\n        H 4 a5bd 3 a5angle 2 a5dihedral\n        C 4 a6bd 3 a6angle 2 a6dihedral\n        H 6 a7bd 4 a7angle 3 a7dihedral\n        C 6 a8bd 4 a8angle 3 a8dihedral\n        H 8 a9bd 6 a9angle 4 a9dihedral\n        C 2 a10bd 3 a10angle 4 a10dihedral\n        H 10 a11bd 2 a11angle 1 a11dihedral\n        H 3 a12bd 2 a12angle 1 a12dihedral\n        a2bd = 1.0853\n        a3bd = 1.3800\n        a4bd = 1.3787\n        a5bd = 1.0816\n        a6bd = 1.3795\n        a7bd = 1.0844\n        a8bd = 1.3872\n        a9bd = 1.0847\n        a10bd = 1.3812\n        a11bd = 1.0838\n        a12bd = 1.0818\n        a3angle = 119.8001\n        a4angle = 121.2010\n        a5angle = 121.0431\n        a6angle = 121.7007\n        a7angle = 122.4276\n        a8angle = 116.6791\n        a9angle = 118.3291\n        a10angle = 117.8148\n        a11angle = 117.3613\n        a12angle = 117.2942\n        a4dihedral = -180\n        a5dihedral = -180\n        a6dihedral = 0\n        a7dihedral = 180\n        a8dihedral = 0\n        a9dihedral = 180\n        a10dihedral = 0\n        a11dihedral = 0\n        a12dihedral = 0\n        '
```
{% endtab %}
{% endtabs %}
