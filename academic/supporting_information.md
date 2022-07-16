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
