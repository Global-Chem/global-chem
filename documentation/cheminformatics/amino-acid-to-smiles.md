# Amino Acid to SMILES



The need for bi-directional for amino acid to SMILES converter seems to be popular in linking data.  To do this we need to borrow code from `CocktailShaker` which marks site points within a string:

```
NCC(NC([*:1])C(NC([*:2])C(O)=O)=O)=O
```

`[*:1]` and `[*:2]` are denoted as site points on the peptide backbone string. We can then install amino acid fragments within these site points.&#x20;

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Convert the Amino Acids Sequence to SMILES**

{% tabs %}
{% tab title="Code" %}
```
# Convert to Amino Acid Sequence to SMILES and Back Again

amino_acid_test = 'RSTEFGHIKLADPQ'
smiles = cheminformatics.amino_acids_to_smiles(amino_acid_test)
print ("SMILES: %s " % smiles)
```
{% endtab %}

{% tab title="Output" %}
```
SMILES: NC(CCCCNC(N)=N)C(NC(CO)C(NC(C(C)([H])O)C(NC(CCC(O)=O)C(NC(CC1=CC=CC=C1)C(NC([H])C(NC(CC1=CNC=N1)C(NC(C(CC)([H])C)C(NC(CCCCN)C(NC(CC(C)C)C(NC(C)C(NC(CC(O)=O)C(NC(C2CCCN2)C(NC(CCC(N)=O)C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O 
```
{% endtab %}
{% endtabs %}

**Convert the SMILES to Amino Acid Sequence**&#x20;

{% tabs %}
{% tab title="Code" %}
```
cheminformatics.smiles_to_amino_acids(smiles)
```
{% endtab %}

{% tab title="Output" %}

{% endtab %}
{% endtabs %}

**Algorithm**&#x20;

{% tabs %}
{% tab title="Amino Acids to SMILES" %}
```
amino_acids_sequence = {
    "A": "C",
    "R": "CCCCNC(N)=N",
    "N": "CC(N)=O",
    "D": "CC(O)=O",
    "B": "CC(O)=O",
    "C": "CS",
    "E": "CCC(O)=O",
    "Q": "CCC(N)=O",
    "Z": "CCC(N)=O",
    "G": "[H]",
    "H": "CC1=CNC=N1",
    "I": "C(CC)([H])C",
    "L": "CC(C)C",
    "K": "CCCCN",
    "M": "CCSC",
    "F": "CC1=CC=CC=C1",
    "P": "C2CCCN2",
    "S": "CO",
    "T": "C(C)([H])O",
    "W": "CCC1=CNC2=C1C=CC=C2",
    "Y": "CC1=CC=C(O)C=C1",
    "V": "C(C)C"
}
```
{% endtab %}

{% tab title="SMILES to Amino Acids" %}
```
amino_acids_sequence = {
            "C" :"A",
            "CCCCNC(N)=N":"R",
            "CC(N)=O":"N",
            "CC(O)=O":"D",
            "CS": "C",
            "CCC(O)=O":"E",
            "CCC(N)=O":"Q",
            "[H]":"G",
            "CC1=CNC=N1" :"H",
            "C(CC)([H])C" :"I",
            "CC(C)C" :"L",
            "CCCCN" :"K",
            "CCSC" :"M",
            "CC1=CC=CC=C1" :"F",
            "C2CCCN2" :"P",
            "CO" :"S",
            "C(C)([H])O" :"T",
            "CCC1=CNC2=C1C=CC=C2" :"W",
            "CC1=CC=C(O)C=C1":"Y",
            "C(C)C" :"V",
}
```
{% endtab %}
{% endtabs %}



So when the user passes in a string of Amino Acids, we build the peptide string in accordance to the length of the string and fill the slots appropriately with amino acid fragments.&#x20;

```
RSTEFGHIKLADPQ
NC(CCCCNC(N)=N)C(NC(CO)C(NC(C(C)([H])O)C(NC(CCC(O)=O)C(NC(CC1=CC=CC=C1)C(NC([H])C(NC(CC1=CNC=N1)C(NC(C(CC)([H])C)C(NC(CCCCN)C(NC(CC(C)C)C(NC(C)C(NC(CC(O)=O)C(NC(C2CCCN2)C(NC(CCC(N)=O)C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O
```

To convert the SMILES string back into an amino acid sequence a reverse list is created and a regex pattern is used to search for the slots within a SMILES string.&#x20;

```
pattern = re.compile('NC\(.*?\)C\(', flags=re.MULTILINE)
```

Then replaces the fragment with the amino acid.&#x20;

```
NC(CCCCNC(N)=N)C(NC(CO)C(NC(C(C)([H])O)C(NC(CCC(O)=O)C(NC(CC1=CC=CC=C1)C(NC([H])C(NC(CC1=CNC=N1)C(NC(C(CC)([H])C)C(NC(CCCCN)C(NC(CC(C)C)C(NC(C)C(NC(CC(O)=O)C(NC(C2CCCN2)C(NC(CCC(N)=O)C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O
RSTEFGHIKLADPQ
```

**References**

1. Sharif, Suliman. “Cocktail Shaker: An Open Source Drug Expansion and Enumeration Library for Peptides.” _Journal of Open Source Software_, vol. 5, no. 52, Aug. 2020, p. 1992. _joss.theoj.org_, https://doi.org/10.21105/joss.01992.
