# GlobalChem

**GlobalChem Internal Network**

class `GlobalChem` was designed where each `Node` is stored as a flat list in the attribute of the class like so as well as itself, `GlobalChem` gets installed as a `dummy` node:

```
class GlobalChem(object):
   
   __NODES__ = [
      'global_chem': Node,
      'emerging_perfluoroalkyls': EmergingPerFluoroAlkyls,
      'montmorillonite_adsorption': MontmorilloniteAdsorption,
   ]
```

If we look at the `initiate`\_`network()`in more detail:

{% tabs %}
{% tab title="Code" %}
```
self.network = {}

self.root_node = root_node

self.network[root_node] = {
  "node_value": Node(
        root_node,
        self.__NODES__[root_node].get_smiles(),
        self.__NODES__[root_node].get_smarts()
   ),
   "children": [],
   "parents": [],
   "name": root_node
 }
```
{% endtab %}

{% tab title="Network Diagram" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 7.53.56 PM.png>)
{% endtab %}
{% endtabs %}

The network gets added with both `parents` and `children` lists. This will act as a linked list of sorts where each node can have as many parents and as many children. Press "Network Diagram" to see the red circle as the node and the purple as the parent and child lists.

Choosing this decision was based on flexibility. There are rules that accommodate this type of software design.&#x20;

**Rule 1**

```
There must be one root node. 
```

For now to start off the graph networks we need to define a center. Somewhere to start and then place all the other nodes. That being said we don't start off we multiple nodes, this was just to make the data structure easier to implement but could be subject to start from several places.

The code can be initialized like so:

```
gc = GlobalChem()
gc.initiate_network()
```

Adding a node is the root of the algorithm and defined here as well as the source code but first there is another rule that must be met:

**Rule 2**

```
There must be a connection from one node to the next. 
Non-Connected Nodes in space cannot exist.
```

For example there cannot be a situation like this:

![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.09.04 PM.png>)

```
Algorithm

1.) First Find the Node in the Network
2.) Add the Child to the Parent
3.) Add the Child to the Network
4.) Add the Parent to the Child
5.) Add the Parent to the Network
```

{% tabs %}
{% tab title="First Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.24.03 PM (2).png>)
{% endtab %}

{% tab title="Second Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.26.06 PM.png>)
{% endtab %}

{% tab title="Third Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.24.51 PM (2).png>)
{% endtab %}

{% tab title="Fourth Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.27.53 PM.png>)
{% endtab %}

{% tab title="Fifth Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.24.51 PM (1).png>)
{% endtab %}
{% endtabs %}

Users can construct their own graphs from here and whatever they fit. This is what makes the lookup so fast instead of traversing through directories since it's only serving as a key:value pair.&#x20;

&#x20;The algorithm is pretty efficient in keep tracking but there is no order which is the trick - another part of the algorithm will handle that which can be called with the function:

```
gc.build_global_chem_network()
```

To build the network first global makes a each python file that is a correspond object into a reversed list:

```
['emerging_perfluoroalkyls', 'environment','global_chem']
```

Then the parent object is popped off:

```
['emerging_perfluoroalkyls', 'environment'] | global_chem
```

Then proceeds into the `node_function` where the next item which is environment is added as a child into `global_chem`. This continues until the length of thew list reaches a value of 1 which is the leaf node.&#x20;

**Deep Graph Layer Network**

The deep layer network follows similar design paradigms but is adjusted a little bit. There must be a root node of 1 which marks as your "input" node. This revisits the case of a similar rule we established before. The fundamental rule is&#x20;

```
Algorithm:
    1.) Fetch all the parents of the current layer
    2.) Add all the node children to the parents.
    3.) Increase the deep layer count.
    4.) Add all the parents to the children.
```

{% tabs %}
{% tab title="First Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 8.59.25 PM.png>)
{% endtab %}

{% tab title="Second Step & Third Step" %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 9.08.33 PM.png>)
{% endtab %}

{% tab title="Fourth Step " %}
![](<../.gitbook/assets/Screen Shot 2022-02-23 at 9.09.33 PM.png>)
{% endtab %}
{% endtabs %}

