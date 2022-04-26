# Sunbursting a SMILES List



{% tabs %}
{% tab title="First  Step" %}
```
for row in self.compounds:

    smiles = row['smiles']
    molecule = Chem.MolFromSmiles(smiles)

    for index, chemical_list in enumerate(self.patterns):

        matches = []

        for name, pattern in chemical_list.items():

            smarts_mol = Chem.MolFromSmarts(pattern)

            try:
                substructs = molecule.GetSubstructMatches(smarts_mol)
                if substructs:
                    matches.append(name)
            except:
                pass

row['%s:matches' % self.record_labels[index]] = matches
```
{% endtab %}

{% tab title="Second Step" %}
```
self.answer_layers = [ [], [] ]

for row in self.compounds:

    for index, chemical_list in enumerate(self.patterns):

        label = self.record_labels[index]
        functional_groups = row['%s:matches' % label]

        if len(functional_groups) > 0:
            for j in functional_groups:
                self.answer_layers[0].append(label)
                self.answer_layers[1].append(j)
```
{% endtab %}
{% endtabs %}

**Third Step**

Plotly uses a series of labels, values, and parents to determine that has to be passed into the sunburst. To do that, we set the root of the parent to be empty and then initial record labels:

The total count of the center of the sunburst needs to be how many were the total sum of the matches of functional groups per record in the set. Then the match numbers of the record to the set.&#x20;

```
first_layer_labels = ['Machine 1'] + self.record_labels
first_layer_parents = [""] + ['Machine 1'] * len(self.record_labels)
first_layer_values = [sum(paper_values.values())] + list(paper_values.values())
```

The second layer is less tricky where you only have to count the number of functional groups that matched within that set.

```
second_layer_labels = [i for i in self.total_values.keys()]
second_layer_parents = [i.split('-')[0].strip() for i in self.total_values.keys()]
second_layer_values = [i for i in self.total_values.values()]
```

Flatten to one list:

```
self.labels = first_layer_labels + second_layer_labels
self.parents = first_layer_parents + second_layer_parents
self.values = first_layer_values + second_layer_values
```

Plot it:

{% tabs %}
{% tab title="Code" %}
```
fig =go.Figure()

        fig.add_trace(go.Sunburst(
            labels=self.labels,
            parents=self.parents,
            values=self.values,
            insidetextorientation='radial',
            #     marker={"colors": colors_charge},
            domain=dict(column=0),
            opacity=0.9,
        ))

        fig.update_layout(
            margin = dict(t=5, l=5, r=5, b=5),
            height=600,
            width=600,
            grid= dict(columns=1, rows=1)
        )
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/newplot - 2022-02-24T120845.390 (5).png>)
{% endtab %}
{% endtabs %}
