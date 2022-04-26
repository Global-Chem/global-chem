# Convert to Networkx Object



`Networkx` is a common graphing tool in python used for graph objects. It is only natural that we also have an extension to convert into networkx objects. This improves interoperability of our software to be included any pipelines. Folk can build chemical graph networks in our program then move into networkx for more processing if they wish.&#x20;

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Convert GlobalChem Graph to NetworkX Graph**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)
network = gc.network
networkx_graph = cheminformatics.convert_to_networkx(network)
print (networkx_graph.nodes.data())
```
{% endtab %}

{% tab title="Output" %}
```
[('global_chem', {}), ('', {}), ('medicinal_chemistry', {}), ('cannabinoids', {}), ('warheads', {}), ('electrophillic_warheads_for_kinases', {}), ('common_warheads_covalent_inhibitors', {}), ('rings', {}), ('phase_2_hetereocyclic_rings', {}), ('rings_in_drugs', {}), ('iupac_blue_book_rings', {}), ('scaffolds', {}), ('iupac_blue_book', {}), ('privileged_scaffolds', {}), ('common_r_group_replacements', {}), ('narcotics', {}), ('pihkal', {}), ('schedule_four', {}), ('schedule_one', {}), ('schedule_two', {}), ('schedule_three', {}), ('schedule_five', {}), ('organic_synthesis', {}), ('bidendate_phosphine_ligands', {}), ('nickel_ligands', {}), ('protecting_groups', {}), ('amino_acid_protecting_groups', {}), ('solvents', {}), ('common_organic_solvents', {}), ('proteins', {}), ('kinases', {}), ('braf', {}), ('braf_inhibitors', {}), ('privileged_kinase_inhibitors', {}), ('formulation', {}), ('excipients', {}), ('biopharmaceutics_class_three', {}), ('cimetidine_acyclovir', {}), ('materials', {}), ('polymers', {}), ('common_monomer_repeating_units', {}), ('clay', {}), ('montmorillonite_adsorption', {}), ('environment', {}), ('emerging_perfluoroalkyls', {}), ('interstellar_space', {}), ('miscellaneous', {}), ('amino_acids', {}), ('regex_patterns', {}), ('vitamins', {}), ('open_smiles', {})]
```
{% endtab %}
{% endtabs %}

**Algorithm**&#x20;

For the conversion, the algorithm is pretty simple:

```
1.) Establish all the nodes
2.) Add all child connections
3.) Add all parent connections
```

We initialize the object:

```
self.networkx_graph = nx.Graph()
```

{% tabs %}
{% tab title="First Step" %}
```
for node_key, node_value in self.network.items():
    self.networkx_graph.add_node(node_value['name'])
```
{% endtab %}

{% tab title="Second Step" %}
```
for node_key, node_value in self.network.items():

   parent_name = node_value['name']
   children = node_value['children']

   for child in children:
       self.networkx_graph.add_edge(child, parent_name)
```
{% endtab %}

{% tab title="Third Step" %}
```
for node_key, node_value in self.network.items():
    parents = node_value['parents']
    child = node_value['name']

    for parent in parents:
        self.networkx_graph.add_edge(child, parent)
```
{% endtab %}
{% endtabs %}

And then return the graph! As simple as that.&#x20;
