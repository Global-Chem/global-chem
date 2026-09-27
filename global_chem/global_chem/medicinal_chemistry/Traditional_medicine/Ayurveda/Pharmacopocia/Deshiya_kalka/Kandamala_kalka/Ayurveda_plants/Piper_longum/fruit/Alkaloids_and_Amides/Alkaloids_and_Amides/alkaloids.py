# Alkaloids and Amides Data for Piper longum
# Chemical structures represented using SMILES, SMARTS patterns, and PubChem IDs

# 1. SMILES Dictionary (Compound Name -> SMILES)
alkaloids_and_amides_smiles = {
    'Isobutyramide': 'CC(C)C(=O)N',
    'Pellitorine': 'CCCCC/C=C/C=C/C(=O)NCC(C)C',
    'Piperlonguminine': 'CC(C)CNC(=O)/C=C/C=C/C1=CC2=C(C=C1)OCO2',
    'Piperine': 'C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3',
    'Piperlongumine': 'COC1=CC(=CC(=C1OC)OC)/C=C/C(=O)N2CCC=CC2=O',
    'Aristolodione': 'CN1C2=CC3=CC=CC=C3C4=C2C(=CC(=C4OC)O)C(=O)C1=O',
    'Norcepharadione B': 'COC1=C(C2=C3C(=C1)C(=O)C(=O)NC3=CC4=CC=CC=C42)OC',
    'Pipernonaline': 'C1CCN(CC1)C(=O)/C=C/CCCC/C=C/C2=CC3=C(C=C2)OCO3',
    'Piperundecalidine': 'C1CCN(CC1)C(=O)/C=C/C=C/CCCC/C=C/C2=CC3=C(C=C2)OCO3',
    'Guineensine': 'CC(C)CNC(=O)/C=C/C=C/CCCCCC/C=C/C1=CC2=C(C=C1)OCO2',
}

# 2. Key Substructures & Functional Groups SMARTS Patterns
alkaloids_and_amides_smarts_patterns = {
    'Amide_Core': '[CX3](=O)[NX3]',
    'Methylenedioxy_Group': 'c1ccc2c(c1)OCO2',
    'Piperidine_Ring': 'C1CCNCC1',
    'Isobutylamide_Moiety': 'CC(C)CNC(=O)',
}

# 3. Complete Dataset Representation (Dictionary of Dictionaries)
alkaloids_and_amides_data = {
    'Isobutyramide': {
        'smiles': 'CC(C)C(=O)N',
        'smarts': '[CX3](=O)[NX3]',
        'pubchem_id': 68424,
    },
    'Pellitorine': {
        'smiles': 'CCCCC/C=C/C=C/C(=O)NCC(C)C',
        'smarts': 'CCCCC/C=C/C=C/C(=O)NCC(C)C',
        'pubchem_id': 5318516,
    },
    'Piperlonguminine': {
        'smiles': 'CC(C)CNC(=O)/C=C/C=C/C1=CC2=C(C=C1)OCO2',
        'smarts': 'CC(C)CNC(=O)/C=C/C=C/c1cc2c(cc1)OCO2',
        'pubchem_id': 5320621,
    },
    'Piperine': {
        'smiles': 'C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3',
        'smarts': 'C1CCN(CC1)C(=O)/C=C/C=C/c2cc3c(cc2)OCO3',
        'pubchem_id': 638024,
    },
    'Piperlongumine': {
        'smiles': 'COC1=CC(=CC(=C1OC)OC)/C=C/C(=O)N2CCC=CC2=O',
        'smarts': 'COc1cc(cc(c1OC)OC)/C=C/C(=O)N2CCC=CC2=O',
        'pubchem_id': 637858,
    },
    'Aristolodione': {
        'smiles': 'CN1C2=CC3=CC=CC=C3C4=C2C(=CC(=C4OC)O)C(=O)C1=O',
        'smarts': 'CN1c2cc3ccccc3c4c2c(cc(c4OC)O)C(=O)C1=O',
        'pubchem_id': 184116,
    },
    'Norcepharadione B': {
        'smiles': 'COC1=C(C2=C3C(=C1)C(=O)C(=O)NC3=CC4=CC=CC=C42)OC',
        'smarts': 'COc1c(c2c3c(=c1)C(=O)C(=O)NC3=cc4ccccc42)OC',
        'pubchem_id': 189168,
    },
    'Pipernonaline': {
        'smiles': 'C1CCN(CC1)C(=O)/C=C/CCCC/C=C/C2=CC3=C(C=C2)OCO3',
        'smarts': 'C1CCN(CC1)C(=O)/C=C/CCCC/C=C/c2cc3c(cc2)OCO3',
        'pubchem_id': 9974595,
    },
    'Piperundecalidine': {
        'smiles': 'C1CCN(CC1)C(=O)/C=C/C=C/CCCC/C=C/C2=CC3=C(C=C2)OCO3',
        'smarts': 'C1CCN(CC1)C(=O)/C=C/C=C/CCCC/C=C/c2cc3c(cc2)OCO3',
        'pubchem_id': 44453654,
    },
    'Guineensine': {
        'smiles': 'CC(C)CNC(=O)/C=C/C=C/CCCCCC/C=C/C1=CC2=C(C=C1)OCO2',
        'smarts': 'CC(C)CNC(=O)/C=C/C=C/CCCCCC/C=C/c2cc3c(cc2)OCO3',
        'pubchem_id': 6442405,
    },
}