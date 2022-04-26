# GlobalChem DNA



Having a useful connector into DNA would be be of an advantage especially an easy conversion to SMILES or SMARTS for comparing potential relations between each other.&#x20;

**Import the Package**

```
from global_chem_extensions import GlobalChemExtensions

bioinformatics = GlobalChemExtensions().bioinformatics()
```

**Initialize the DNA Molecule**&#x20;

{% tabs %}
{% tab title="Code" %}
```
dna_molecule bioinformatics.intialize_globalchem_dna(
    dna_sequence='AAAATGTGTGTG'
)
```
{% endtab %}

{% tab title="Output" %}
```
<global_chem_extensions.entities.dna.dna.GlobalChemDNA at 0x7f81547f5550>
```
{% endtab %}
{% endtabs %}

&#x20;**Convert to SMILES**

{% tabs %}
{% tab title="Code" %}
```
dna_molecule.convert_to_smiles()
```
{% endtab %}

{% tab title="Output" %}
```
Nc1nc2nc(C3CC(OP(=O)([O-])OCC4OC(Cc5c[nH]c(=O)[nH]c5=O)CC4OP(=O)([O-])OCC4OC(c5nc6nc(N)[nH]c(=O)c6[nH]5)CC4OP(=O)([O-])OCC4OC(Cc5c[nH]c(=O)[nH]c5=O)CC4OP(=O)([O-])OCC4OC(c5nc6nc(N)[nH]c(=O)c6[nH]5)CC4OP(=O)([O-])OCC4OC(Cc5c[nH]c(=O)[nH]c5=O)CC4OP(=O)([O-])OCC4OC(c5nc6nc(N)[nH]c(=O)c6[nH]5)CC4OP(=O)([O-])OCC4OC(Cc5c[nH]c(=O)[nH]c5=O)CC4OP(=O)([O-])OCC4OC(c5nc6ncnc(N)c6[nH]5)CC4OP(=O)([O-])OCC4OC(c5nc6ncnc(N)c6[nH]5)CC4OP(=O)([O-])OCC4OC(c5nc6ncnc(N)c6[nH]5)CC4OP(=O)([O-])OCC4OC(c5nc6ncnc(N)c6[nH]5)CC4OP(=O)([O-])[O-])C(COP(=O)([O-])OC4CCOC4CO)O3)[nH]c2c(=O)[nH]1
```
{% endtab %}
{% endtabs %}

**Convert to SMARTS**

{% tabs %}
{% tab title="Code" %}
```
dna_molecule.convert_to_smarts()
```
{% endtab %}

{% tab title="Output" %}
```
[#7]-[#6]1:[#7]:[#6]2:[#7]:[#6](-[#6]3-[#6]-[#6](-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]-[#6]5:[#6]:[#7H]:[#6](=[#8]):[#7H]:[#6]:5=[#8])-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6](-[#7]):[#7H]:[#6](=[#8]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]-[#6]5:[#6]:[#7H]:[#6](=[#8]):[#7H]:[#6]:5=[#8])-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6](-[#7]):[#7H]:[#6](=[#8]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]-[#6]5:[#6]:[#7H]:[#6](=[#8]):[#7H]:[#6]:5=[#8])-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6](-[#7]):[#7H]:[#6](=[#8]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]-[#6]5:[#6]:[#7H]:[#6](=[#8]):[#7H]:[#6]:5=[#8])-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6]:[#7]:[#6](-[#7]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6]:[#7]:[#6](-[#7]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6]:[#7]:[#6](-[#7]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]-[#6]4-[#8]-[#6](-[#6]5:[#7]:[#6]6:[#7]:[#6]:[#7]:[#6](-[#7]):[#6]:6:[#7H]:5)-[#6]-[#6]-4-[#8]-[#15](=[#8])(-[#8-])-[#8-])-[#6](-[#6]-[#8]-[#15](=[#8])(-[#8-])-[#8]-[#6]4-[#6]-[#6]-[#8]-[#6]-4-[#6]-[#8])-[#8]-3):[#7H]:[#6]:2:[#6](=[#8]):[#7H]:1
```
{% endtab %}
{% endtabs %}

**Visualize the DNA Strand**

{% tabs %}
{% tab title="Code" %}
```
dna_molecule.visualize_strand()
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/Screen Shot 2022-04-07 at 11.01.39 PM.png>)
{% endtab %}
{% endtabs %}

**Visualize and Label the DNA Strand**

{% tabs %}
{% tab title="Code" %}
```
dna_molecule.label_feature(
    start=0,
    end=3,
    label='Start'
)

dna_molecule.label_feature(
    start=4,
    end=8,
    label='Middle'
)


dna_molecule.label_feature(
    start=4,
    end=8,
    label='Middle 2',
    color='Blue'
)

dna_molecule.label_feature(
    start=8,
    end=12,
    label='End',
    color='Green'
)

dna_molecule.visualize_strand()
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/Screen Shot 2022-04-07 at 11.04.39 PM.png>)
{% endtab %}
{% endtabs %}

**Save the DNA Image to a file**

```
dna.save_to_image() 
```
