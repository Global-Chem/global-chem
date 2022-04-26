# Principal Component Analysis



Principal component analysis is a well studied technique in identifying feature characteristics on a set of data and the variance of the data. The more some data aligns with a certain range we can determine as a cluster. In the context of chemistry, when we have a large set of `SMILES` data if we would like to identify major features that collect groups of molecules that might not be as obvious as before.&#x20;

**Imports**

```
from global_chem import GlobalChem
from global_chem_extensions import GlobalChemExtensions

gc = GlobalChem()
cheminformatics = GlobalChemExtensions().cheminformatics()
```

**Conduct PCA Analysis**

{% tabs %}
{% tab title="Code" %}
```
gc = GlobalChem()
gc.build_global_chem_network(print_output=False, debugger=False)
smiles_list = list(gc.get_node_smiles('schedule_one').values())

mol_ids = cheminformatics.node_pca_analysis(
            smiles_list,
            morgan_radius = 1,
            bit_representation = 512,
            number_of_clusters = 5,
            number_of_components = 0.95,
            random_state = 0,
            save_file=False,
            return_mol_ids=True,
)
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/hellooo (1).gif>)
{% endtab %}
{% endtabs %}

**Algorithm**&#x20;

To accomplish this, we need numbers and not a list of strings. We apply the the conversion of SMILES to fingerprints. First we need to capture a local chemical environment and convert to a series of numbers. Morgan Fingerprinting was used initially with something called the "Morgan Radius" which is a series of discrete integers that represent how many bonds to iterate to look.

![](<../.gitbook/assets/Screen Shot 2022-02-24 at 8.56.46 AM.png>)

The morgan radius at 0 only looks at the direct connections to the atoms and as you increase the radius you look the atoms connections that are connected to your direct connections. Doing this we can start to evaluate the chemical environment.&#x20;

```
arr = np.zeros((0,))

fp = AllChem.GetMorganFingerprintAsBitVect(
    mol, 
    self.morgan_radius, 
    nBits=self.bit_representation
)

DataStructs.ConvertToNumpyArray(fp, arr)
```

The conversion into the fingerprints means how much feature space you want to capture. Usually a bit length has a defined amount and the more numbers than in theory the more different features. Although, there is a level of redundancy and too wild of feature space. The default for ours is 512 which enough to get an overall look at a list of SMILES data without knowing anything in the first place.

For Principal Component Analysis, we want to evaluate the variance of the chemical data and start grouping clusters together that correspond to similar bit patterns. Do this we apply a vector to a series of data and a component orthogonal to the original vector.&#x20;

&#x20;

![](<../.gitbook/assets/Screen Shot 2022-02-24 at 9.12.27 AM.png>)

We then evaluate the variance of the data with respect to the vector. Now how you place your vectors is where the variability comes in the data. And how many vectors you want to add is also subjective. The resultant plot is a linear combination of the variance of the data which reduces the data into two main categories: outliers and normal values.&#x20;

To add some context in the terms of bit representation:

```
molecules_list = [Chem.MolFromSmiles(i) for i in self.smiles_list]
fingerprints_list = np.array([self.mol2fparr(m) for m in molecules_list])

pca = PCA(n_components=self.number_of_components)
chemicalspace = pca.fit_transform(fingerprints_list)
```

But this is not enough because the chemical space from the pca analysis on large sets of SMILES can have a lot of variance. To further gain some insight from the pca we employ a k-means clustering on the subspace:

```
kmean = KMeans(n_clusters=self.number_of_clusters, random_state=self.random_state)
kmean.fit(fingerprints_list)
kmeanc = [colormaps[i] for i in kmean.labels_]
```

To get a hierarchal view of the major core scaffolds in our ligand set. The result produced can produce variety of clusters depending on the data:

![](<../.gitbook/assets/Screen Shot 2022-02-24 at 9.23.28 AM.png>)

**References**

1. Jolliffe, I. T., editor. “Principal Component Analysis and Factor Analysis.” _Principal Component Analysis_, Springer, 2002, pp. 150–66. _Springer Link_, https://doi.org/10.1007/0-387-22440-8\_7.
2. Ding, Chris, and Xiaofeng He. “_K_ -Means Clustering via Principal Component Analysis.” _Twenty-First International Conference on Machine Learning  - ICML ’04_, ACM Press, 2004, p. 29. _DOI.org (Crossref)_, https://doi.org/10.1145/1015330.1015408.CloseDeleteEdit
3.  Morgan, H. L. “The Generation of a Unique Machine Description for Chemical Structures-A Technique Developed at Chemical Abstracts Service.” Journal of Chemical Documentation, vol. 5, no. 2, May 1965, pp. 107–13. DOI.org (Crossref), https://doi.org/10.1021/c160017a018**.**

