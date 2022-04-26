# One-Hot Encoding SMILES

One-Hot encoding is the processing of moving something that is built on some sort of language can be indexed appropriately with numbers. In machine learning, it is common to convert strings to numeric characters so the machine can work with the numbers. Most often there is an encoder to convert the strings, perform the machine learning, and then decode the numbers to check the performance of the machine learning.

{% embed url="https://towardsdatascience.com/basic-molecular-representation-for-machine-learning-b6be52e9ff76" %}

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Convert A List of SMILES to GlobalChem Encoder**

{% tabs %}
{% tab title="Code" %}
```
smiles_list = list(gc.get_node_smiles('pihkal').values())

encoded_smiles = cheminformatics.encode_smiles(smiles_list, max_length=200)

print ('Encoded SMILES: %s' % encoded_smiles[0])

decoded_smiles = cheminformatics.decode_smiles(encoded_smiles)

print ('Decoded SMILES: %s' % decoded_smiles[0])
```
{% endtab %}

{% tab title="Output" %}
```
Encoded SMILES: [[0. 0. 0. ... 0. 0. 0.]
 [0. 0. 0. ... 0. 0. 0.]
 [0. 0. 0. ... 0. 0. 0.]
 ...
 [0. 0. 0. ... 0. 0. 0.]
 [0. 0. 0. ... 0. 0. 0.]
 [0. 0. 0. ... 0. 0. 0.]]
Decoded SMILES: CCC(N)CC1=CC(=C(OC)C(=C1)OC)OC
```
{% endtab %}
{% endtabs %}

**Algorithm**&#x20;

In the code there is an index mapping of the the SMILES to the index and then convert to numpy array. Currently, the mapping is set to this to include polymeric SMILES, atom classes, and virtual atoms.&#x20;

```
__SMILES_MAPPING__ = [
        ' ',
        '#', '%', '(', ')', '+', '-', '.', '/',
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        '=', '@',
        'A', 'B', 'C', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P',
        'R', 'S', 'T', 'V', 'X', 'Z',
        '[', '\\', ']',
        'a', 'b', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'p', 'r', 's',
        't', 'u',
        '&', ':', '*'
    ]
```
