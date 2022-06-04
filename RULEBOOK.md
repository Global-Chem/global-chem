
GlobalChem Rulebook
===================

Hi and welcome to the official Rulebook of GlobalChem to how the data is prepared and curated. To maintain data integreity for a common language of chemistry
we need to adhere to rules that make it easier for the human to pronounce the words with maximial information with least amount of letters or numbers. 

### IUPAC

#### Rule 1: Remove Numeric Characters in favour of Preferred Natural Name Except if Number is Natural. 

In majority of the cases we want to remove numbers or go in favour of names that don't have numbers embeddded with the IUPAC key. That's because we want something
more natural sounding like. 

For example, let's take a simple example of cocaine:

```python
3-Tropanylbenzoate-2-carboxylic acid methyl ester
```

which is more commonly known as `cocaine`, so we prefer the general language of the common person.

In some cases where we report narcotics such as `2C` class family which their chemicals are most commonly known as their class '2c' followed by the one-letter derivative to represent the analogue. 

```python
2CI, 2CE
```

#### Rule 2: Remove Stereochemistry (Trans and Cis) from both SMILES and IUPAC which can be determined algorithmically. 

We want to remove any track of trans or cis since that can be determined algorithmically in the smiles and then installed in the IUPAC string. This reduces
down the information and makes the language more approachable. Let's say we evaluate common cannabinoid compounds where the trans is reported. 

```python
delta-9-(trans)-tetrahydrocannabinol
``

We remove the `trans`:

```python
delta-9-tetrahydrocannabinol
```

The SMILES also follows a similar protocol:

```python
CCCCCc1cc(c2c(c1)OC([C@H]3[C@H]2C=C(CC3)C)(C)C)O
```

Where the stereochemistry is reported with `@` or `@@` symbol to represent `R` and `S` configurations. We opt to remove the square brackets surrounding
the hydrogen and the hydrogen to:

```python
CCCCCc1cc(c2c(c1)OC(C3C2C=C(CC3)C)(C)C)O
```

