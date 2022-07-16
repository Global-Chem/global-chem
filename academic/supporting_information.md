## Data Collection

References and associated compound lists are selected based on the interests of the scientific contributors.  This should include consideration of relevance to the scientific community. To authenticate and validate SMILES strings we employ interoperability tools to other cheminformatic software to verify it's usability.

The IUPAC/SMILES strings may be abstracted in a variety of methods:

- For IUPAC naming we opted for naming things as they were reported in the literature. If no names were available, then we opted to find a natural name to fill the slot.

- For IUPAC naming, we chose to reduce the complexity of the name by opting to remove as much stereochemistry as made sense. 

- For Polymer IUPAC, the site points were omitted from the name and some of the nomenclature adjusted for preferred names
over traditional. For example: 'yl' to mark site points for polymer connections was removed in favor of reduced english complexity. Site points are marked with a virtual atom that can be installed into the SMILES string with the character '*'.

-  For simple molecules one representation of the SMILES can be directly translated using visual 
inspection. This is typically appropriate for compounds at the beginning of a reported list that contain the most common denominator rings. 

- For complex molecules the image can be redrawn in the free version of ChemDraw and then translated into SMILES. 

- For sources where the SMILES are written and the IUPAC is not known the SMILES are translated into ChemDraw and the name retrieved. 
Note that some of the names may be modified based on human inspection in favor of preferred names. 

- In the case of radicals, some SMILES were adjusted to remove the radical chemical feature as they serve as connection points. However in some cases the radical component was maintained, especially in the case of IUPAC blue book common substituents or instellar space where radicals are more unknown and not as well explored.

- SMARTS strings were adapted from the SMILES using RDKit [@Landrum:2019-5]. 

- Stereochemistry which is represented as '@' and '@@' as R and S respectively are removed from the SMILES in order to reduce complexity.


## Cheminformatics Test

Global-Chem parsed through seven different tools with majority being successful minus diamond represented with an '&' [Reference Here] and fails with all software including RDKit except the GlobalChem Encoder does account for it. (it is not clear) The percentage of passing is as follows: RDKit 100% [Reference Here], DeepSMILES 99.25% [Reference Here], PartialSMILES 85.7% [Reference Here] , SELFIES 100% [Reference Here], MolVS 98.5% [Reference Here], PySMILES 99.8% [Reference Here]. PartialSMILES proved to be the most robust acceptance/rejection tool in identifying misrepresentations of SMILES. 

| Software        | Number of Failed Compounds | Example of Failed SMILES                                   |
|-----------------|----------------------------|------------------------------------------------------------|
| RDKit           | 0                          |                                                            |
| SELFIES         | 0                          |                                                            |
| Indigo          | 8                          | 'CC([Si](C1=CC=CC=C1)C2=CC=CC=C2)(C)C'                     |
| PySMILES        | 5                          | '[a].[Na+].[K+].[Mg+2].[Ca+2].[Cl-]'                       |
| DeepSMILES      | 8                          | 'OC(C1=CC=CC=C1N2)C2=O', 'C1(N2C(CN=C3)=NN=C2)=C3C=CC=C1'  |    
| MolVS           | 24                         | 'n1ccnc1', 'HF', 'O=N1CCCCC1'                              |
| PartialSMILES   | 337                        | '[CH]C', '[N]=[N+]=[N-]'                                   |

<p align="center">
  <i>Table 2: Intereoperability Results of Common Moleculesfrom Global-Chem against different cheminformatic software</i>
</p>

Indigo's encoder was pretty robust and their software allows for a lot of inteoperabiltiy with different software tools (i.e pdf data parsing of SMILES), when faced with the tert-butyldiphenylsilyl protecting group and the SMILES string with the `Si` is not wrapped in a square brackets that specify an element that doesn't have a complete valence shell. For PySMILES, the inclusion of the '[a]' denoting aromaticity for an "aromatic salt" in the database couldn't be processed. Some other encoders have encoded for an aromaticity keyword as specified in the Daylight Technical Documentation [Reference Here].  DeepSMILES was interesting because it failed on specific functional groups as shown in the example with an oxindole and triazolodiazepines that had complex small branch complexities and moetieies that it didn't foresee existing. MolVS had some interesting results where imidazole (and derivatives) failed probably because it expected for a hydrogen perhaps to be explicity stated due to it's varying protonation states. Hydrofluoric acid was something I was expecting but again the hydrogen actually needed to be enforced with a [H] which is not as intuitive. PartialSMILES proved to be the most robust eluding to SMILES that were partially complete and rejected by their criteria. Failures included a ethyl radical and a azido complex stemming from the interstellarspace molecules. 
