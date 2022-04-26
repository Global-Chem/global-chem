# Psi4Parser & Orbital Visualizer



Psi4 is rising as a popular open source QM tool and it makes sense for there to be a connection to the output and potentially visualizing cube files for molecule orbitals.&#x20;

First initialize the quantum chemistry package

```
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChemExtensions().quantum_chemistry()
```

**Initialize the Parser**

```
parser = gc.get_psi4_parser()
```

**Read in a Psi4 Out File**

```
parser = parser(input_file='test.out')
```

**Extract Energy Values**

The parser get's main energetic terms from the psi4 out file while excluding other things that are probably not needed.

{% tabs %}
{% tab title="Code" %}
```
parser.get_energy_contributions()
```
{% endtab %}

{% tab title="Output" %}
```
{
   'nuclear_repulsion_energy': 9.353841051132276, 
   'one_electron_energy': -123.40409212008015, 
   'two_electron_energy': 37.99248795110696, 
   'total_energy': -76.05776311784092
}
```
{% endtab %}
{% endtabs %}

**Extract Dipole Moments**

Extract the Dipole Moments in their relative x, y, and z-axis.

{% tabs %}
{% tab title="Code" %}
```
parser.get_dipole_contributions()
```
{% endtab %}

{% tab title="Output" %}
```
{
   'nuclear_dipole_moment': ['0.0000', '0.0000', '0.9422'], 
   'electronic_dipole_moment': ['0.0000', '0.0000', '-0.1602'],
  'dipole_moment_atomic_units': ['0.0000', '0.0000', '0.7821', '0.7821'], 
   'dipole_moment_debye': ['0.0000', '0.0000', '1.9878', '1.9878']
}
```
{% endtab %}
{% endtabs %}

**Extract Frequencies**

If the frequencies or the harmonics were run then the parser is able to handle extracting all the modes both real and not. Psi4 usually only shows the IR active modes and truncates the rest.&#x20;

{% tabs %}
{% tab title="Code" %}
```
parser.get_frequencies()
```
{% endtab %}

{% tab title="Output" %}
```
['0.0000i', '0.0000i', '0.0000', '0.0000', '0.0000', '0.0001', '1752.9915', '4127.2470', '4227.2071']
```
{% endtab %}
{% endtabs %}

**Visualize Cube Files**&#x20;

Visualize the cube prop files using moly to get a look at the HOMO/LUMO of your molecule

{% tabs %}
{% tab title="Code" %}
```
parser.moly_plot_molecular_orbital(
    'Psi_a_1_1-A1.cube'
)
```
{% endtab %}
{% endtabs %}
