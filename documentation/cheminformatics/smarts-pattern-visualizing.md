# SMARTS Pattern Visualizing

This mainly comes from the REST service provided by the SMARTS plus group at the University of Hamburg:

{% embed url="https://smarts.plus" %}

**Imports**

```
from global_chem_extensions import GlobalChemExtensions

cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Visualize the SMILES**

```
molecule = cheminformatics.visualize_smarts(
'[#6]-[#7]1-[#6](-[#6]=[#6]-[#6]-1=[#8])=[#8]'
)

molecule
```

**Algorithm**&#x20;

The code is pretty simple using the base url and passing in the SMARTs string:

```
 self.base_url = "https://smarts.plus/smartsview/download_rest?"
 self.url = self.base_url + f"smarts={self.smarts_pattern}"
```

To run the code:

```
from global_chem_extensions import GlobalChemExtensions

molecule = GlobalChemExtensions.visualize_smarts(
'[#6]-[#7]1-[#6](-[#6]=[#6]-[#6]-1=[#8])=[#8]'
)

molecule
```

From here we can analyze different parts of the SMARTs string and test different things out.&#x20;

![](<../.gitbook/assets/Screen Shot 2022-03-31 at 8.06.03 AM (2).png>)

The breakdown of the string:
