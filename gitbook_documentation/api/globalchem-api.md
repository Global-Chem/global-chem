# GlobalChem API



First initialize the class must be initialized:

```
gc = GlobalChem()
```

**Print the Network**

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

**Check the available nodes in GlobalChem:**

{% tabs %}
{% tab title="Code" %}
```
nodes_list = gc.check_available_nodes()
print (nodes_list)
```
{% endtab %}

{% tab title="Output" %}
```
'emerging_perfluoro_alkyls', 
'montmorillonite_adsorption', 
'common_monomer_repeating_units',
 'electrophilic_warheads_for_kinases',
```
{% endtab %}
{% endtabs %}

**Build the GlobalChem Network and Print it Out**

{% tabs %}
{% tab title="Code" %}
```
gc.build_global_chem_network(
    print_output=True,          # Print the network out
    debugger=False,             # For Developers mostly to see all node values
)
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
},
```
{% endtab %}
{% endtabs %}

The algorithm uses a series of parents/children to connect nodes instead of "edges" as in traditional graph networks. This just makes it easier to code if the graph database lives as a 1-dimensional with lists of parents and children's connected in this fashion.

**Fetch a Node**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network()
node = gc.get_node('emerging_perfluoroalkyls')
print (node)
```
{% endtab %}

{% tab title="Output" %}
```
{
    'node_value': <global_chem.global_chem.Node object at 0x117fee210>, 
    'children': [], 
    'parents': ['emerging_perfluoroalkyls'], 
    'name': 'emerging_perfluoroalkyls'
}
```
{% endtab %}
{% endtabs %}

**Fetch the IUPAC:SMILES/SMARTS Data from the Node**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network()
smiles = gc.get_node_smiles('emerging_perfluoroalkyls')
smarts = gc.get_node_smarts('emerging_perfluoroalkyls')

print ("Length of Perfluoroalkyls: %s " % len(smiles))
```
{% endtab %}

{% tab title="Output" %}
```
Length of Perfluoroalkyls: 27 
```
{% endtab %}
{% endtabs %}

**Fetch All Data From the Network**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
print(gc.get_all_smiles())
print(gc.get_all_smarts())
print(gc.get_all_names())
```
{% endtab %}

{% tab title="Output" %}
```
C(=O)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)O
[#6](=[#8])(-[#6](-[#6](-[#6](-[#6](-[#6](-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])-[#8]
perfluorohexanoic acid
```
{% endtab %}
{% endtabs %}

**Remove a Node from the Network**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)
gc.remove_node('emerging_perfluoroalkyls')
print (list(gc.network.keys()))
```
{% endtab %}

{% tab title="Output" %}
```
['global_chem', '', 'organic_synthesis', 'solvents', 
'common_organic_solvents', 'protecting_groups',
'amino_acid_protecting_groups', 'materials', 'polymers',
'common_monomer_repeating_units', 'clay',
'montmorillonite_adsorption', 'proteins', 'kinases',
'scaffolds', 'privileged_kinase_inhibitors', 
'braf', 'braf_inhibitors', 'miscellaneous', 
'vitamins', 'open_smiles', 'amino_acids', 
'regex_patterns', 'environment', 'narcotics', 'schedule_one',
'schedule_four', 'schedule_five', 'pihkal', 'schedule_two', 
'schedule_three', 'interstellar_space', 'medicinal_chemistry', 
'cannabinoids', 'privileged_scaffolds', 'iupac_blue_book',
'common_r_group_replacements', 'warheads',
'electrophillic_warheads_for_kinases',
'common_warheads_covalent_inhibitors', 
'rings',
'phase_2_hetereocyclic_rings', 
'iupac_blue_book_rings', 
'rings_in_drugs'
 ]

```
{% endtab %}
{% endtabs %}

**Fetch a SMILES By the IUPAC**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
definition = gc.get_smiles_by_iupac(
    'benzene',                          
    return_network_path=False,          # Return the last found network path
    return_all_network_paths=True,      # Return all the found network paths 
)
print (definition)
```
{% endtab %}

{% tab title="Output" %}
```
{'C1=CC=CC=C1': 'global_chem.miscellaneous.open_smiles'}
```
{% endtab %}
{% endtabs %}

**Set & Get Node value in the Network**

Node values are here to act as a way for the user to elect a list or a subset list perhaps in heural learning or it could really be anything the user wants however they construct their network.&#x20;

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network(print_output=True, debugger=False)
gc.set_node_value('emerging_perfluoroalkyls', {'some_data': ['bunny']})
print (gc.get_node_value('emerging_perfluoroalkyls'))
```
{% endtab %}

{% tab title="Output" %}
```
{'some_data': ['bunny']}
```
{% endtab %}
{% endtabs %}

**To Create Your Own Chemical Graph Network (GN) And Check the Values**

{% tabs %}
{% tab title="Code" %}
```
from global_chem import GlobalChem

gc = GlobalChem(verbose=False)
gc.initiate_network()
gc.add_node('global_chem', 'common_monomer_repeating_units')
gc.add_node('common_monomer_repeating_units','electrophilic_warheads_for_kinases')
values = gc.get_node_smiles('common_monomer_repeating_units')

print (values)
```
{% endtab %}

{% tab title="Output" %}
```
'3′-bromo-2-chloro[1,1′:4′,1′′-terphenyl]-4,4′′':
 'ClC1=CC=CC=C1C2=CC=C(C3=CC=CC=C3)C(Br)=C2'

```
{% endtab %}
{% endtabs %}

**Creating a Deep Layer Chemical Graph Networks (DGN) & Print it Out:**

This is for the more advanced users of building networks and how to manage sets of layers.&#x20;

{% tabs %}
{% tab title="Code" %}
```
# Create a Deep Layer Network

gc = GlobalChem()
gc.initiate_deep_layer_network()
gc.add_deep_layer(
    [
        'emerging_perfluoroalkyls',
        'montmorillonite_adsorption',
        'common_monomer_repeating_units'
    ]
)
gc.add_deep_layer(
    [
        'common_warheads_covalent_inhibitors',
        'privileged_scaffolds',
        'iupac_blue_book'
    ]
)

gc.print_deep_network()
```
{% endtab %}

{% tab title="Output" %}
```
                                      ┌common_warhead_covalent_inhibitors
            ┌emerging_perfluoroalkyls─├privileged_scaffolds
            │                         └iupac_blue_book
            │                           ┌common_warhead_covalent_inhibitors
global_chem─├montmorillonite_adsorption─├privileged_scaffolds
            │                           └iupac_blue_book
            │                               ┌common_warhead_covalent_inhibitors
            └common_monomer_repeating_units─├privileged_scaffolds
                                            └iupac_blue_book

```
{% endtab %}
{% endtabs %}

**Compute a Common Score**

Common Score Algorithm:

1. Datamine the current graph network of GlobalChem
2. Get the object weights of each mention&#x20;
3. Determine the mention weight
4. Sum the Weight's and that is how common the molecule is.&#x20;

The higher the value the higher the common score tied with it's IUPAC name.&#x20;

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)
gc.compute_common_score('benzene', verbose=True)
```
{% endtab %}

{% tab title="Output" %}
```
GlobalChem Common Score: 7.139921294271971
```
{% endtab %}
{% endtabs %}

**To TSV**

&#x20;The network returned in all CSV format for interoperability for web application development mostly but can also be used to search.&#x20;

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.to_tsv('global_chem.tsv')
```
{% endtab %}
{% endtabs %}

