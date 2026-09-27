# Lignans Data for Piper longum
# Chemical structures represented using SMILES, SMARTS patterns, and PubChem IDs

# 1. SMILES Dictionary
lignans_smiles = {
    'Pluviatilol': 'COC1=C(C=CC(=C1)[C@@H]2[C@@H]3CO[C@H]([C@@H]3CO2)C4=CC5=C(C=C4)OCO5)O',
    'Fargesin': 'COC1=C(C=C(C=C1)[C@H]2[C@H]3CO[C@@H]([C@H]3CO2)C4=CC5=C(C=C4)OCO5)OC',
}

# 2. Key Substructures & Functional Groups SMARTS Patterns
lignans_smarts_patterns = {
    'Furofuran_Core': 'C1OCC2C1COC2',
    'Methylenedioxy_Group': 'c1ccc2c(c1)OCO2',
    'Methoxyphenyl_Group': 'COc1ccccc1',
}

# 3. Complete Dataset Representation
lignans_data = {
    'Pluviatilol': {
        'smiles': 'COC1=C(C=CC(=C1)[C@@H]2[C@@H]3CO[C@H]([C@@H]3CO2)C4=CC5=C(C=C4)OCO5)O',
        'smarts': 'COc1c(ccc(c1)[C@@H]2[C@@H]3CO[C@H]([C@@H]3CO2)c4cc5c(cc4)OCO5)O',
        'pubchem_id': 70695727,
    },
    'Fargesin': {
        'smiles': 'COC1=C(C=C(C=C1)[C@H]2[C@H]3CO[C@@H]([C@H]3CO2)C4=CC5=C(C=C4)OCO5)OC',
        'smarts': 'COc1c(cc(cc1)[C@H]2[C@H]3CO[C@@H]([C@H]3CO2)c4cc5c(cc4)OCO5)OC',
        'pubchem_id': 10926754,
    },
}