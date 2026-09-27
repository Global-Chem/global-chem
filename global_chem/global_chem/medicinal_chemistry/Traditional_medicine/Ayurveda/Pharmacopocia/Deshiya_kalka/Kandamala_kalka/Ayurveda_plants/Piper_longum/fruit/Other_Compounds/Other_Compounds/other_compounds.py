# Other Compounds Data for Piper longum
# Chemical structures represented using SMILES, SMARTS patterns, and PubChem IDs

# 1. SMILES Dictionary
other_compounds_smiles = {
    '1-Methylhexyl acetate': 'CCCCCC(C)OC(=O)C',
}

# 2. Key Substructures & Functional Groups SMARTS Patterns
other_compounds_smarts_patterns = {
    'Ester_Group': '[CX3](=O)[OX2H0][CX4]',
    'Acetate_Moiety': 'CC(=O)O',
}

# 3. Complete Dataset Representation
other_compounds_data = {
    '1-Methylhexyl acetate': {
        'smiles': 'CCCCCC(C)OC(=O)C',
        'smarts': 'CCCCCC(C)OC(=O)C',
        'pubchem_id': 80018,
    },
}