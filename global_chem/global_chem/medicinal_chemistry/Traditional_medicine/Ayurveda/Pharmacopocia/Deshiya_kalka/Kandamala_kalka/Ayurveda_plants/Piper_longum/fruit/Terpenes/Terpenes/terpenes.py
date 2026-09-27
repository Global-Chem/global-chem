# Terpenes Data for Piper longum
# Chemical structures represented using SMILES, SMARTS patterns, and PubChem IDs

# 1. SMILES Dictionary
terpenes_smiles = {
    'Alpha-Terpineol': 'CC1=CCC(CC1)C(C)(C)O',
    'Alpha-pinene': 'CC1=CCC2CC1C2(C)C',
    'Limonene': 'CC1=CCC(CC1)C(=C)C',
    'Linalool': 'CC(=CCCC(C)(C=C)O)C',
    'Myrcene': 'CC(=CCCC(=C)C=C)C',
    'Sabinene': 'CC(C)C12CCC(=C)C1C2',
    'Terpinolene': 'CC1=CCC(=C(C)C)CC1',
    'Beta-caryophyllene': 'C/C/1=C\\CCC(=C)[C@H]2CC([C@@H]2CC1)(C)C',
    'Cadinane': 'C[C@H]1CC[C@H]2[C@H](CC[C@H]([C@@H]2C1)C(C)C)C',
    'Humulene': 'C/C/1=C\\CC(/C=C/C/C(=C/CC1)/C)(C)C',
    'Zingiberene': 'CC1=CC[C@@H](C=C1)[C@@H](C)CCC=C(C)C',
    'alpha-Muurolene': 'CC1=C[C@@H]2[C@H](CC1)C(=CC[C@H]2C(C)C)C',
}

# 2. Key Substructures & Functional Groups SMARTS Patterns
terpenes_smarts_patterns = {
    'Isoprene_Unit': 'CC(=C)C=C',
    'Cyclohexene_Ring': 'C1CCC=CC1',
    'Bicyclo_Heptane_System': 'C1C2CCC1C2',
}

# 3. Complete Dataset Representation
terpenes_data = {
    'Alpha-Terpineol': {
        'smiles': 'CC1=CCC(CC1)C(C)(C)O',
        'smarts': 'CC1=CCC(CC1)C(C)(C)O',
        'pubchem_id': 17100,
    },
    'Alpha-pinene': {
        'smiles': 'CC1=CCC2CC1C2(C)C',
        'smarts': 'CC1=CCC2CC1C2(C)C',
        'pubchem_id': 6654,
    },
    'Limonene': {
        'smiles': 'CC1=CCC(CC1)C(=C)C',
        'smarts': 'CC1=CCC(CC1)C(=C)C',
        'pubchem_id': 22311,
    },
    'Linalool': {
        'smiles': 'CC(=CCCC(C)(C=C)O)C',
        'smarts': 'CC(=CCCC(C)(C=C)O)C',
        'pubchem_id': 6549,
    },
    'Myrcene': {
        'smiles': 'CC(=CCCC(=C)C=C)C',
        'smarts': 'CC(=CCCC(=C)C=C)C',
        'pubchem_id': 31253,
    },
    'Sabinene': {
        'smiles': 'CC(C)C12CCC(=C)C1C2',
        'smarts': 'CC(C)C12CCC(=C)C1C2',
        'pubchem_id': 18818,
    },
    'Terpinolene': {
        'smiles': 'CC1=CCC(=C(C)C)CC1',
        'smarts': 'CC1=CCC(=C(C)C)CC1',
        'pubchem_id': 11463,
    },
    'Beta-caryophyllene': {
        'smiles': 'C/C/1=C\\CCC(=C)[C@H]2CC([C@@H]2CC1)(C)C',
        'smarts': 'C/C/1=C\\CCC(=C)[C@H]2CC([C@@H]2CC1)(C)C',
        'pubchem_id': 5281515,
    },
    'Cadinane': {
        'smiles': 'C[C@H]1CC[C@H]2[C@H](CC[C@H]([C@@H]2C1)C(C)C)C',
        'smarts': 'C[C@H]1CC[C@H]2[C@H](CC[C@H]([C@@H]2C1)C(C)C)C',
        'pubchem_id': 9548708,
    },
    'Humulene': {
        'smiles': 'C/C/1=C\\CC(/C=C/C/C(=C/CC1)/C)(C)C',
        'smarts': 'C/C/1=C\\CC(/C=C/C/C(=C/CC1)/C)(C)C',
        'pubchem_id': 5281520,
    },
    'Zingiberene': {
        'smiles': 'CC1=CC[C@@H](C=C1)[C@@H](C)CCC=C(C)C',
        'smarts': 'CC1=CC[C@@H](C=C1)[C@@H](C)CCC=C(C)C',
        'pubchem_id': 92776,
    },
    'alpha-Muurolene': {
        'smiles': 'CC1=C[C@@H]2[C@H](CC1)C(=CC[C@H]2C(C)C)C',
        'smarts': 'CC1=C[C@@H]2[C@H](CC1)C(=CC[C@H]2C(C)C)C',
        'pubchem_id': 12306047,
    },
}