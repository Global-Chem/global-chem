# Deep Layer Scatter

The Deep Layer Scatter Plot is still in prototype phase but really good for seeing where the package is going in terms of lexical network graphs. We implement the Parallel Coordinates from Plotly:

{% embed url="https://plotly.com/python/parallel-coordinates-plot" %}

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Scatter your Graph**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.initiate_deep_layer_network()
gc.add_deep_layer(
    [
        'emerging_perfluoroalkyls',
        'montmorillonite_adsorption',
        'common_monomer_repeating_units',
        'vitamins',
        'schedule_two',
    ]
)
gc.add_deep_layer(
    [
        'common_warheads_covalent_inhibitors',
        'privileged_scaffolds',
        'iupac_blue_book',
        'schedule_four',
        'schedule_one',
    ]
)

gc.print_deep_network()e_file=True)
```
{% endtab %}

{% tab title="Output" %}
```
                                      ┌common_warheads_covalent_inhibitors
                                      ├privileged_scaffolds
            ┌emerging_perfluoroalkyls─├iupac_blue_book
            │                         ├schedule_four
            │                         └schedule_one
            │                           ┌common_warheads_covalent_inhibitors
            │                           ├privileged_scaffolds
            ├montmorillonite_adsorption─├iupac_blue_book
            │                           ├schedule_four
            │                           └schedule_one
            │                               ┌common_warheads_covalent_inhibitors
            │                               ├privileged_scaffolds
global_chem─├common_monomer_repeating_units─├iupac_blue_book
            │                               ├schedule_four
            │                               └schedule_one
            │         ┌common_warheads_covalent_inhibitors
            │         ├privileged_scaffolds
            ├vitamins─├iupac_blue_book
            │         ├schedule_four
            │         └schedule_one
            │             ┌common_warheads_covalent_inhibitors
            │             ├privileged_scaffolds
            └schedule_two─├iupac_blue_book
                          ├schedule_four
                          └schedule_one
```
{% endtab %}
{% endtabs %}

**Algorithm**

Here we want to visualize deep layer networks and add colour, weights of lines etc for a more sophisticated data visualization.

In the algorithm we first determine the max depth of the network object. While keep tracking of how many nodes are in each deep layer. For example, the `Network Graph` we have **3** layers&#x20;

{% tabs %}
{% tab title="Code" %}
```
layer = 1
node_counts = []

for node_key, node_value in self.deep_layer_network.items():

    if int(node_value['layer']) > layer:
        layer = int(node_value['layer'])

for i in range(1, layer + 1):

    counts = 0

    for node_key, node_value in self.deep_layer_network.items():

        if int(node_value['layer']) == i:
            counts += 1

    node_counts.append(counts)
```
{% endtab %}

{% tab title="Network Graph" %}
![](<../.gitbook/assets/Screen Shot 2022-02-24 at 11.20.01 AM.png>)
{% endtab %}
{% endtabs %}

Now we have to convert our GlobalChem Network object into a Plotly Style to get it ready for `Parcats` plotting.&#x20;

```
Algorithm:
    1.) Establish the First Layer as Root & Establish Second Layer Concomitantly.
    2.) Determine Lexical Symmetry or Anti-Symmetry & Add Continued Layers
    3.) Prepare Plotly Python Data Structure
```

In the third step, the dimensions refers to the a list that will house the list of dictionary objects. Each `category` is a layer, where it's values correspond to a combination of lexical keys. For example let's say we create a three node graph with 3 layers:&#x20;

```
Layer 1: ['global_chem']
Layer 2: ['global_chem - emerging_perfluroalkyls' ]
Layer 3: ['global_chem - emerging_perfluroalkyls - schedule_one' ]
```

{% tabs %}
{% tab title="First Step" %}
```
if i == 1:
    first_layer = [ nodes[0]['name'] ] * reduce((lambda x, y: x * y), self.node_counts)
    layer_object.append(first_layer)

    values = nodes[0]['children']

    for value in values:
        second_layer += [ value ] * self.node_counts[1]

    layer_object.append(second_layer)
```
{% endtab %}

{% tab title="Second Step" %}
```
for node in nodes:

    if not symmetric:
        layer += [ node['name'] ]

    else:
        layer += [ node['name'] ]  * self.node_counts[i - 1]

if not symmetric:

    layer = layer  * self.node_counts[i - 1]

layer_object.append(layer)
```
{% endtab %}

{% tab title="Third Step" %}
```
dimensions = []
colors = []

for index, layer in enumerate(layers_labels):

    categories = {}

    categories['label'] = layer
    categories['values'] = layer_object[index]

    if index == 1:
        colors += [ 'purple' ] * self.node_counts[0]

    elif index == 2:
        colors += [ 'lightseagreen' ] * self.node_counts[1]

    else:
        colors += [ 'gold' ] * self.node_counts[index - 1]


    categories['categoryorder'] = 'category ascending'
```
{% endtab %}
{% endtabs %}

{% hint style="info" %}
Lexical Symmetry refers to the connection order, where you add layers in a symmetric and anti-symmetric fashion to obtain all combinations
{% endhint %}

With our dimensions we then pass into `Parcats`

{% tabs %}
{% tab title="Code" %}
```
fig = go.Figure(go.Parcats(
            dimensions=dimensions,
            arrangement='freeform',
            line = dict(
                color = colors,
            ),
            hoveron='color',
            hoverinfo='count+probability',
            labelfont={'size': 20, 'family': 'Times'},
            tickfont={'size': 15, 'family': 'Times'},
        ))

        fig.update_layout(
            title_text="Deep Layer Scatter",
            title_font=dict(size=12, family='Arial'),
            template='simple_white',
            xaxis_tickformat = 'i',
            bargap=0.4, # gap between bars of adjacent location coordinates,
            height=800,
            width=1700
        )
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/newplot - 2022-02-24T113615.213.png>)
{% endtab %}
{% endtabs %}
