# Validation



To perform validation of the SMILES, we'll be passing it through a series of software to determine what fails and what doesn't.&#x20;

The full performance of the GlobalChem against these different lists can be found below:



| Software       | SMILES Passed |
| -------------- | ------------- |
| RDKit          | 100%          |
| SELFIES        | 100%          |
| PySMILES       | 99.8%         |
| Partial SMILES | 85.7%         |
| DeepSMILES     | 99.25%        |
| MolVS          | 98.50%        |

Users can also validate their own SMILES across different software:

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Validate a list SMILES across different software**

{% tabs %}
{% tab title="Code" %}
```
gc.build_global_chem_network()
smiles_list = list(gc.get_node_smiles('emerging_perfluoroalkyls').values())

sucesses, failures = cheminformatics.verify_smiles(
    smiles_list,
    rdkit=True, 
    partial_smiles=True,
    return_failures=True,
    pysmiles=True,
    molvs=True
)

total = len(sucesses) + len(failures)
print ("Percantage of Accepted SMILES: %s" % ((len(sucesses) / total) * 100))
```
{% endtab %}
{% endtabs %}
