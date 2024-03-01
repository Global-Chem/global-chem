# Welcome to GlobalChem Contributing Guide! 

Thank you for investing your time in contributing to our project! Any contribution you make will be reflected on either the node contribution list or the extension contribution list. In this guide you will get an overview of how to get started in contributing into GlobalChem.

## Contributing to GlobalChem

### Filing an Issue

GlobalChem is a collection of lists and most likely if you published an academic paper highlighting structural trends for small molecules or gathered data from some other source you would
most likely want people to use your data so you can create a node and add it to our network. 

First file an issue, place your initials of username or name followed by an issue number that is in respect to you. For me it is:

```
SS-1, SS-2, SS-3
```

But for let's say one of other developers it is:

```
BL-1, Bl-2, BL-3
```

File an issue and reference the paper or data source you would like to add. Let's use the Common Excipients as our example: https://doi.org/10.1002/jps.24643.

### Adding your Node

Let's say we want to add common excipients for these biopharmaceutics class III of compounds of Cimetidine and acyclovir. We add data in two forms IUPAC/SMILES, and IUPAC/SMARTS. Each node is represented as a class
object like so:

We can write a file called the node like so and each node is written in the exact same way.

```python
class CimetidineAndAcyclovir(object):

    def __init__(self):

        self.name = 'cimetidine_and_acyclovir'
```

The class is camel case where the name attribute is separated by a `_`. So in this paper there is a list of chemical functional groups relevant to paper. What we do
is read the functional groups and translate them down into a language we can write as english, IUPAC. Now if you don't know how to begin, you can actually
use this free service by ChemDraw online:

https://chemdrawdirect.perkinelmer.cloud/js/sample/index.html#

Redraw the molecules into here and then click this button::

<img width="941" alt="Screen Shot 2022-05-08 at 9 09 51 AM" src="https://user-images.githubusercontent.com/11812946/167297711-f0f434a3-65bd-437d-af0b-17fc9891b7b6.png">

When that happens you will get:

<img width="272" alt="Screen Shot 2022-05-08 at 9 09 41 AM" src="https://user-images.githubusercontent.com/11812946/167297718-d66918fd-c3dd-49c0-bf26-604e34b77ef0.png">

A name for it. Now you can add that name to the list like so:

```python

class CimetidineAndAcyclovir(object):

    def __init__(self):

        self.name = 'cimetidine_and_acyclovir'

    @staticmethod
    def get_smiles():

        smiles = {
            'hexane-1,2,3,4,5,6-hexaol': '',
        }
        
    return
```

Now part of GlobalChem's benefit is that we started to move into more common language natural names rather than IUPAC defined rules. This is something we can speak. Traditionally, I want to remove numbers entirely in the dataset which will be a big curation day. So as we progress further that's something we want to preserve. If you put that name in google, chances are you will find synonyms like a thesauras and what we are looking for is something we can pronounce like `mannitol`. I can say it out loud. Now you can also use ChemDraw to download the SMILES or write it yourself and place it in the string.

```python

class CimetidineAndAcyclovir(object):

    def __init__(self):

        self.name = 'cimetidine_and_acyclovir'

    @staticmethod
    def get_smiles():

        smiles = {
            'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
        }
        
    return
```    

 Do that for the whole list!
 
 #### Bonus
 
If you want to generate the SMARTS string. I usually take the SMILES string and convert to `RDKit` Mol Object and then to SMARTS. If the SMILES is written nicely the SMARTS is written nicely.


```python

from rdkit import Chem

smiles = {
    'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
}

if __name__ == '__main__':

    for k, v in smiles.items():

        print (f"'{k}': '{Chem.MolToSmarts(Chem.MolFromSmiles(v))}',")
```

Commit your file on a branch and we will handle the organization as it gets included into the network!


## Contributing to GlobalChemExtensions

Let's say you have a function or feature you would like to add that you think would be useful. Most often in academia with the wide use of python there are loose scripts floating around that could be of useful functionality but no proper way to distribute it or a medium. If there is a proper distribution most often it is drowned in the noise of publications. To restore some of this useful features, we wet lab scientists are exploring lots of different software to highlight the best and add it to our collection so we can preserve it's functionality. 

### Adding your script

If you have a python file that you would like to add to the `GlobalChemExtensions` then please navigate to the folder structure which looks like this:

<img width="899" alt="Screen Shot 2022-05-08 at 9 20 46 AM" src="https://user-images.githubusercontent.com/11812946/167298165-88d9f636-f689-4237-886b-9eb36874d7e5.png">

Each field is an object with applications that are useful for them:

<img width="949" alt="Screen Shot 2022-05-08 at 9 21 56 AM" src="https://user-images.githubusercontent.com/11812946/167298204-ccc899fa-d739-458c-9de7-7247ce81feb2.png">

The cheminformatics have a plethora of functionality spread from different packages:

<img width="920" alt="Screen Shot 2022-05-08 at 9 22 56 AM" src="https://user-images.githubusercontent.com/11812946/167298232-930dceb4-a3d3-4916-9a3e-f0f3b42820ba.png">

Add your script here and explain your functionality. Submit a PR and we will do a code review. We hope to maintain good quality but also recognize the expedience of academia and will help refactor your code. Not afraid of doing that.

Hope to see more of you guys and meet people across different places as we grow our graph.

Cheers.

