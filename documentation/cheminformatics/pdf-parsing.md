# PDF Parsing



PDF Parsing is going to be handling by a separate package called MolPDF. It has pretty much a fixed template and is used for handling data distribution quickly so we can get a general feel of a molecule list.&#x20;

{% embed url="https://github.com/Sulstice/molpdf" %}

The philosophy is very simple. We don't really need a template but a simple PDF document to transfer data quickly especially in the case of a supplementary information. We use the meta data stored in the PDF to store the SMILES in accordance with the images. That way it makes it easier to mine the data quickly into a python object.&#x20;

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**SMILES to PDF**

{% tabs %}
{% tab title="Code" %}
```
smiles_list = list(gc.get_node_smiles('pihkal').values())

cheminformatics.smiles_to_pdf(
    smiles=smiles_list,
    labels = [],
    file_name = 'molecules.pdf',
    include_failed_smiles = True,
    title = 'MY MOLECULES',
)
```
{% endtab %}

{% tab title="Output" %}
```
Method: 'generate' Time: 4.98 seconds
```
{% endtab %}
{% endtabs %}

**PDF to SMILES**

{% tabs %}
{% tab title="Code" %}
```
molecules = cheminformatics.pdf_to_smiles(
    'molecules.pdf',
)

print (len(molecules))
```
{% endtab %}
{% endtabs %}
