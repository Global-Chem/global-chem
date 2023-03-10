
from rdkit import Chem
from rdkit.Chem import AllChem


aa_key = {
    'Arg': 'R',
    'His': 'H',
    'Lys': 'K',
    'Asp': 'D',
    'Glu': 'E',
    'Ser': 'S',
    'Thr': 'T',
    'Asn': 'N',
    'Gln': 'Q',
    'Cys': 'C',
    'Gly': 'G',
    'Pro': 'P',
    'Ala': 'A',
    'Val': 'V',
    'Ile': 'I',
    'Leu': 'L',
    'Met': 'M',
    'Phe': 'F',
    'Tyr': 'Y',
    'Trp': 'W'
}

smiles = {
    'heptane': 'CCCCCCC',
    '2-methyl heptane': 'CCCCCC(C)C',
    '3-methyl heptane': 'CCCCC(C)CC',
    'octane': 'CCCCCCCC',
    'ethyl cyclohexane': 'CCC1CCCCC1',
    '4-methyl octane': 'CCCCC(C)CCC',
    '2-methyl octane': 'CCCCCCC(C)C',
    '2,5-dimethyl heptane': 'CCC(C)CCC(C)C',
    '3-methyl octane': 'CCCCCC(C)CC',
    '2,2,4-trimethyl heptane': 'CCCC(C)CC(C)(C)C',
    'butyl cyclopentane': 'CCCCC1CCCC1',
    'nonane': 'CCCCCCCCC',
    '3-ethyl-2,5-dimethyl hexane': 'CCC(CC(C)C)C(C)C',
    '1,3-dimethyl benzene': 'CC1=CC(=CC=C1)C',
    'propyl cyclohexane': 'CCCC1CCCCC1',
    '2,5,6-trimethyl octane': 'CCC(C)C(C)CCC(C)C',
    '2-methyl nonane': 'CCCCCCCC(C)C',
    '3-methyl nonane': 'CCCCCCC(C)CC',
    '2,2,4,6,6-pentamethyl heptane': 'CC(CC(C)(C)C)CC(C)(C)C',
    '2,3-dimethyl nonane': 'CCCCCCC(C)C(C)C',
    'decane': 'CCCCCCCCCC',
    '1-methyl-2-propyl cyclohexane': 'CCCC1CCCCC1C',
    '5-ethyl-2,2,3-trimethyl heptane': 'CCC(CC)CC(C)C(C)(C)C',
    '1-ethyl-2-methyl benzene': 'CCC1=CC=CC=C1C',
    '3-methyl decane': 'CCCCCCCC(C)CC',
    '1,2,3-trimethyl benzene': 'CC1=C(C(=CC=C1)C)C',
    '3,7-dimethyl decane': 'CCCC(C)CCCC(C)CC',
    '4-methyl undecane': 'CCCCCCCC(C)CCC',
    'undecane': 'CCCCCCCCCCC',
    '1-ethenyl-2-methyl benzene': 'C=CC1=CC=CC=C1C',
    '1,2-diethyl benzene': 'CCC1=CC=CC=C1CC',
    '1-methyl-4-(2-propenyl) benzene': 'CC1=CC=C(C=C1)C=C(C)C',
    '1-ethyl-2,4-dimethyl benzene': 'CCC1=C(C=C(C=C1)C)C',
    '1-methyl-2-(4-methylpentyl)cyclopentane': 'CC1CCCC1CCCC(C)C',
    '3-methyl undecane': 'CCCCCCCCC(C)CC',
    '1,2,3,5-tetramethyl benzene': 'CC1=CC(=C(C(=C1)C)C)C',
    '2-methyl 1,1*-bicyclohexyl': 'CC1CCCCC1C2CCCCC2',
    '1,6-tridecadiene': 'CCCCCCC=CCCCC=C',
    'dodecane': 'CCCCCCCCCCCC',
    '1-ethenyl-4-ethyl benzene': 'CCC1=CC=C(C=C1)C=C',
    '1-methyl-2-(2-propenyl) benzene': 'CC=CC1=CC=CC=C1C',
    '2,6-dimethyl undecane': 'CCCCCC(C)CCCC(C)C',
    '1,2,3,4-tetrahydro naphthalene': 'C1CCC2=CC=CC=C2C1',
    '2-ethenyl-1,3,5-trimethyl benzene': 'C=CC1=C(C)C=C(C)C=C1C',
    'hexyl cyclohexane': 'CCCCCCC1CCCCC1',
    '3-methyl dodecane': 'CCCCCCCCCC(C)CC',
    '1,2,3,4-tetrahydro-2-methyl naphthalene': 'CC1CCC2=CC=CC=C2C1',
    'tridecane': 'CCCCCCCCCCCCC',
    '2,2,4,4,6,8,8-heptamethyl nonane': 'CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C',
    '1,2,3,4-tetrahydro-5-methyl naphthalene': 'CC1=C2CCCCC2=CC=C1',
    '3-methyl tridecane': 'CCCCCCCCCCC(C)CC',
    '2,6,10-trimethyl dodecane': 'CCC(C)CCCC(C)CCCC(C)C',
    '1,2,3,4-tetrahydro-2,7-diemthyl napthalene': 'CC1CCC2=C(C1)C=C(C=C2)C',
    'tetradecane': 'CCCCCCCCCCCCCC',
    '7-methyl pentadecane': 'CCCCCCCCC(C)CCCCCC',
    '3-methyl tetradecane': 'CCCCCCCCCCCC(C)CC',
    'pentadecane': 'CCCCCCCCCCCCCCC',
    'hexadecane': 'CCCCCCCCCCCCCCCC',
    'heptadecane': 'CCCCCCCCCCCCCCCCC',
}

if __name__ == '__main__':
    # sequence = 'Asp Cys Ser Ser Gly Trp Ser Ser Tyr Glu Gly His Cys Tyr Lys Val Phe Lys Gln Ser Lys Thr Trp Ala Asp Ala Glu Ser Phe Cys Thr Lys Gln Val Asn Gly Gly His Leu Val Ser Leu Glu Ser Ser Gly Glu Ala Asp Phe Val Gly Gln Leu Leu Ala Gln Lys Leu Lys Ser Ala Lys Leu His Val Trp Leu Gly Leu Arg Ala Gln Asn Lys Glu Lys Gln Cys Ser Leu Gln Trp Ser Asp Gly Ser Ser Leu Ser Tyr Glu Asn Trp Leu Glu Glu Glu Ser Lys Lys Cys Leu Gly Val His Leu Glu Thr Gly Phe His Lys Trp Glu Asn Phe Tyr Cys Glu Gln Gln Asp Pro Phe Val Cys Glu Ala Asp Cys Pro Ser Asp Trp Ser Ser Tyr Glu Gly His Cys Tyr Lys Pro Phe Asn Glu Pro Lys Asn Trp Ala Asp Ala Glu Asn Phe Cys Thr Lys Gln His Thr Gly Gly His Leu Val Ser Phe Gln Ser Thr Glu Glu Ala Asp Phe Val Val Lys Leu Ala Phe Gln Thr Phe Asp Tyr Gly Leu Phe Trp Phe Gly Leu Ser Lys Leu Trp Asn Gln Cys Asn Trp Gln Trp Ser Asn Ala Ala Met Leu Lys Tyr Thr Asp Trp Ala Glu Glu Ser Tyr Cys Val Tyr Phe Lys Ser Thr Asn Asn Lys Trp Arg Ser Leu Thr Cys Arg Met Leu Ala Asn Phe Val Cys Glu Phe Gln Ala'
    #
    # sequence = sequence.split(' ')
    #
    # aa = [key[i] for i in sequence]
    # print (''.join(aa))
    #
    print (len(smiles))
    for i, v in smiles.items():

        print ("'%s': '%s'," % (i, Chem.MolToSmarts(Chem.MolFromSmiles(v))))

    # for i, v in smiles.items():
    #
    #     molecule = Chem.MolFromSmiles(v)
    #     bit_string = AllChem.GetMorganFingerprintAsBitVect(
    #         molecule,
    #         2,
    #         nBits=512
    #     ).ToBitString()
    #
    #     print ("'%s': '%s'," % (i, bit_string))



