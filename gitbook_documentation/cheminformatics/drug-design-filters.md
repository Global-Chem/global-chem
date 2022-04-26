# Drug Design Filters

Most common are drug filters in reducing a dataset of large SMILES down. These filters can be cumbersome to find, compare, and code. So we extended `GlobalChemExtensions` to cover that.  We use `RDKit` to fetch the parameters of each molecule.&#x20;

**Imports**

```
from global_chem_extensions import GlobalChemExtensions

cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Filter a Drug by one of the Filters**

{% tabs %}
{% tab title="Code" %}
```
gc.build_global_chem_network()
smiles_list = list(gc.get_node_smiles('emerging_perfluoroalkyls').values())

filtered_smiles = cheminformatics.filter_smiles_by_criteria(
    smiles_list,
    lipinski_rule_of_5=True,
    ghose=False,
    veber=False,
    rule_of_3=False,
    reos=False,
    drug_like=False,
    pass_all_filters=False
)

print (filtered_smiles)
```
{% endtab %}

{% tab title="Output" %}
```
[
   'O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 
   'O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
   'O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
   'O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
   'FC(F)(F)C1(F)OC1(F)F', 
   'O=C(O)C(F)(F)C(F)(F)C(F)(F)F',
   'O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 
   'O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
]

```
{% endtab %}
{% endtabs %}

**Algorithm**&#x20;

**Lipinski**

{% tabs %}
{% tab title="Rules" %}
```
Lipinski:
    Moleculer Weight <= 500
    LogP <= 5
    H-Bond Donor Count <= 5
    H-Bond Acceptor Count <= 10
```
{% endtab %}

{% tab title="Implementation" %}
```
if molecular_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10and rotatable_bonds <= 5:
    lipinski = True
    results["Lipinski Rule of 5"] += 1
```
{% endtab %}
{% endtabs %}

**Ghose**

{% tabs %}
{% tab title="Rules" %}
```
Ghose:
    Molecular weight between 160 and 480
    LogP between -0.4 and +5.6
    Atom count between 20 and 70
    Molar refractivity between 40 and 130
```
{% endtab %}

{% tab title="Implementation" %}
```
if molecular_weight >= 160 and molecular_weight <= 480 and logp >= -0.4 and logp <= 5.6 and number_of_atoms >= 20 and number_of_atoms <= 70 and molar_refractivity >= 40 and molar_refractivity <= 130:
    ghose_filter = True
    results["Ghose Filter"] += 1
```
{% endtab %}
{% endtabs %}

**Veber**

{% tabs %}
{% tab title="Rules" %}
```
Veber:
    Rotatable bonds <= 10
    Topological polar surface area <= 140
```
{% endtab %}

{% tab title="Implementation" %}
```
if rotatable_bonds <= 10 and topological_surface_area_mapping <= 140:
    veber_filter = True
    results["Veber Filter"] += 1
```
{% endtab %}
{% endtabs %}

**Rule of 3**

{% tabs %}
{% tab title="Rules" %}
```
Rule of 3:
    Molecular weight <= 300
    LogP <= 3
    H-bond donor <= 3
    H-bond acceptor count <= 3
    Rotatable bond count <= 3
```
{% endtab %}

{% tab title="Implementation" %}
```
if molecular_weight <= 300 and logp <= 3 and h_bond_donor <= 3 and h_bond_acceptors <= 3 and rotatable_bonds <= 3:
    rule_of_3 = True
    results["Rule of 3 Filter"] += 1
```
{% endtab %}
{% endtabs %}

**Drug-Like (QED)**

{% tabs %}
{% tab title="Rules" %}
```
mass < 400
ring count > 0
rotatable bond count < 5
h-bond donor count <= 5
h-bond acceptor count <= 10
logP < 5
```
{% endtab %}

{% tab title="Implementation" %}
```
if molecular_weight < 400 and num_of_rings > 0 and rotatable_bonds < 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and logp < 5:
    drug_like_filter = True
    results["Drug-like Filter"] += 1
```
{% endtab %}
{% endtabs %}

For a more detailed look, please visit this blog for an application on the `WITHDRAWN` database:
