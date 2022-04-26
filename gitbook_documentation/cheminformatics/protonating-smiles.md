# Protonating SMILES

Local chemical environment is of interest to anyone involved in drug design. The protonation state of a molecule or protein can have a lot of affects of the binding affinity or pose of the molecule as it goes it's reaction. Knowing the protonation state prior to actually building the drug will help us design drugs more accurately and possibly breach into covalent binding as well.&#x20;

To implement a design we lean on the open source package `DimorphiteDL` and have it determine and sanity check the SMILES of a string.  The pka is calculated on the fly determined by the user's input for minimum and maximum ph.&#x20;

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Protonate your SMILES**

{% tabs %}
{% tab title="Code" %}
```
protonation_states = cheminformatics.find_protonation_states(
    ['CC(=O)O'],
    min_ph=4.0,
    max_ph=8.1,
    pka_precision=1.0,
    max_variants=128,
    label_states=False,
)

print (protonation_states)
```
{% endtab %}

{% tab title="Output" %}
```
{'CC(=O)O': ['CC(=O)O', 'CC(=O)[O-]']}
```
{% endtab %}
{% endtabs %}

**References**

1.) Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.
