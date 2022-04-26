# Quick Start

Here you will learn how to install the initial package, build the initial network, print it out, and perform some analysis.&#x20;

## Install the Package(s)

**GlobalChem** is the graph network that has no dependencies and it's functionality built within the main object.

**GlobalChemExtensions** has a dependency network that is too cumbersome to deal with but has additional functionality for cheminformaticians (or anyone really) to perform analysis on chemical data including the GlobalChem Graph network.&#x20;

The best way to interact with our API is to use one our official library distributed on PyPi

{% tabs %}
{% tab title="Python" %}
```
# Install via pip
pip install global-chem 

# Install the Extension Package
pip install global-chem[extensions]
```
{% endtab %}
{% endtabs %}

**Additional Dependency Features**

Not everyone wants to install everything into their local environment which can be a very hefty especially something as large as `GlobalChem`To combat this we partitioned some of the applications dependencies into different package dependencies that can be installed with the extra function from `setuptools`. Please refer to the master extension list about which app depends on where.

* **web\_server** - any web applications that depend on the flask app
* **validation** - any SMILES validation dependencies
* **bioinformatics** - any software components that depend on the `biopython` or `dna_resolver`
* **machine\_learning** - any software packages that rely on encoding or learning or preparing data for machine learning.&#x20;
* **pdf** - PDF Generation from SMILES and back gain which require image dependencies.&#x20;
* **plotting -** Any libraries that pertain to plotting mostly in `plotly` or `bokeh`
* **all -** All the extensions

{% tabs %}
{% tab title="Code" %}
```
pip install global-chem[web_server]
pip install global-chem[validation]
pip install global-chem[bioinformatics]
pip install global-chem[machine_learning]
pip install global-chem[pdf]
pip install global-chem[plotting]
pip install global-chem[all]
```
{% endtab %}
{% endtabs %}

{% hint style="info" %}
**Good to know:** global-chem-extensions dependencies are not linked to any specific versions in hopes for flexibility of other development environments.&#x20;
{% endhint %}

## GlobalChem: Building the GlobalChem Graph Network

To build the GlobalChem Graph network we first import the package, initialize the class, and call the function :

{% tabs %}
{% tab title="Code" %}
```
from global_chem import GlobalChem
gc = GlobalChem()
gc.build_global_chem_network(print_output=True)
```
{% endtab %}

{% tab title="Output" %}
```
'global_chem': {
    'children': [
        'environment',
        'miscellaneous',
        'organic_synthesis',
        'medicinal_chemistry',
        'narcotics',
        'interstellar_space',
        'proteins',
        'materials'
    ],
    'name': 'global_chem',
    'node_value': <global_chem.global_chem.Node object at 0x10f60eed0>,
    'parents': []
}, etc.
```
{% endtab %}
{% endtabs %}

To make it pretty:

{% tabs %}
{% tab title="Code" %}
```
gc.print_globalchem_network()
```
{% endtab %}

{% tab title="Output" %}
```

                                ┌solvents─common_organic_solvents
             ┌organic_synthesis─└protecting_groups─amino_acid_protecting_groups
             │          ┌polymers─common_monomer_repeating_units
             ├materials─└clay─montmorillonite_adsorption
             │                            ┌privileged_kinase_inhibtors
             │                            ├privileged_scaffolds
             ├proteins─kinases─┌scaffolds─├iupac_blue_book_substituents
             │                 │          └common_r_group_replacements
             │                 └braf─inhibitors
             │              ┌vitamins
             │              ├open_smiles
             ├miscellaneous─├amino_acids
             │              └regex_patterns
global_chem──├environment─emerging_perfluoroalkyls
             │          ┌schedule_one
             │          ├schedule_four
             │          ├schedule_five
             ├narcotics─├pihkal
             │          ├schedule_two
             │          └schedule_three
             ├interstellar_space
             │                    ┌cannabinoids
             │                    │         ┌electrophillic_warheads_for_kinases
             │                    ├warheads─└common_warheads_covalent_inhibitors
             └medicinal_chemistry─│      ┌phase_2_hetereocyclic_rings
                                  └rings─├iupac_blue_book_rings
                                         └rings_in_drugs
                                         
```
{% endtab %}
{% endtabs %}

## GlobalChem Extensions Access Nodes and Perform a PCA  and Radial Analysis:

Let's have some fun. Let's access a node and perform some PCA Analysis. We want to test whether an object functional groups share some similarity some arbitrary features and try to determine what those features specifically are. This will help understand features of relevance for small molecules.&#x20;

We are going to look at the list of the molecules in pihkal because it's a pretty comprehensive list of what's on the drug market currently published on the wikipedia page. This will help us identify&#x20;

**PCA Analysis**

{% tabs %}
{% tab title="Code" %}
```javascript
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)
smiles_list = list(gc.get_node_smiles('pihkal').values())

GlobalChemExtensions().node_pca_analysis(smiles_list, save_file=False)
```
{% endtab %}

{% tab title="Plot" %}
![](<.gitbook/assets/Screen Shot 2022-02-24 at 8.43.47 AM.png>)

If we look at the PCA analysis in more detail, we can see that the machine has assigned fingerprints of the benzodioxole and the secondary amine with long alkyl chain as primary features for this cluster. If you hover over dots you can fetch other clusters that seem as primary core scaffolds within this dataset.&#x20;
{% endtab %}
{% endtabs %}

**Radial Analysis**

Let's have a look at how a list of emerging perfluoroalkyls to the rest of the nodes in the network using a radial analysis. For more details on the Radial Analysis algorithm please head over to the page.&#x20;

{% tabs %}
{% tab title="Code" %}
```javascript
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)

smiles_list = list(gc.get_node_smiles('pihkal').values())
GlobalChemExtensions().sunburst_chemical_list(smiles_list, save_file=False)
```
{% endtab %}

{% tab title="Plot" %}
![](<.gitbook/assets/newplot - 2022-02-24T085008.217.png>)



If we have a quick look at the list of Pihkal, we can see that they are very similar to the covalent warheads. More nodes, and more in-ferment can be made about the data and up to the user verify :).&#x20;

![](<.gitbook/assets/newplot - 2022-02-20T095100.543.png>)
{% endtab %}
{% endtabs %}

**Enjoy**

Read more of the documentation or just start playing around with the data. This data takes some time to digest so patience is necessary when building you're own networks as well. Happy cheminformatics.&#x20;



****

##
