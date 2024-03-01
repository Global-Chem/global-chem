
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
    'desogestrel': 'CCC12CC(=C)C3C(C1CCC2(C#C)O)CCC4=CCCCC34',
    'dienogest': 'CC12CCC3=C4CCC(=O)C=C4CCC3C1CCC2(CC#N)O',
    'drospirenone': 'CC12CCC(=O)C=C1C3CC3C4C2CCC5(C4C6CC6C57CCC(=O)O7)C',
    'estradiol': 'CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O',
    'estradiol valerate': 'CCCCC(=O)OC1CCC2C1(CCC3C2CCC4=C3C=CC(=C4)O)C',
    'estriol': 'CC12CCC3C(C1CC(C2O)O)CCC4=C3C=CC(=C4)O',
    'estrone': 'CC12CCC3C(C1CCC2=O)CCC4=C3C=CC(=C4)O',
    'ethinyl estradiol': 'CC12CCC3C(C1CCC2(C#C)O)CCC4=C3C=CC(=C4)O',
    'etonogestrel': 'CCC12CC(=C)C3C(C1CCC2(C#C)O)CCC4=CC(=O)CCC34',
    'levonorgestrel': 'CCC12CCC3C(C1CCC2(C#C)O)CCC4=CC(=O)CCC34',
    'mestranol': 'CC12CCC3C(C1CCC2(C#C)O)CCC4=C3C=CC(=C4)OC',
    'norelgestromin': 'CCC12CCC3C(C1CCC2(C#C)O)CCC4=CC(=NO)CCC34',
    'norethisterone': 'CC12CCC3C(C1CCC2(C#C)O)CCC4=CC(=O)CCC34',
    'norethisterone acetate': 'CC(=O)OC1(CCC2C1(CCC3C2CCC4=CC(=O)CCC34)C)C#C',
    'norethynodre': 'CC12CCC3C(C1CCC2(C#C)O)CCC4=C3CCC(=O)C4',
    'norgestimate': 'CCC12CCC3C(C1CCC2(C#C)OC(=O)C)CCC4=CC(=NO)CCC34',
    'progesterone': 'CC(=O)C1CCC2C1(CCC3C2CCC4=CC(=O)CCC34C)C',
}

if __name__ == '__main__':
    # sequence = 'Asp Cys Ser Ser Gly Trp Ser Ser Tyr Glu Gly His Cys Tyr Lys Val Phe Lys Gln Ser Lys Thr Trp Ala Asp Ala Glu Ser Phe Cys Thr Lys Gln Val Asn Gly Gly His Leu Val Ser Leu Glu Ser Ser Gly Glu Ala Asp Phe Val Gly Gln Leu Leu Ala Gln Lys Leu Lys Ser Ala Lys Leu His Val Trp Leu Gly Leu Arg Ala Gln Asn Lys Glu Lys Gln Cys Ser Leu Gln Trp Ser Asp Gly Ser Ser Leu Ser Tyr Glu Asn Trp Leu Glu Glu Glu Ser Lys Lys Cys Leu Gly Val His Leu Glu Thr Gly Phe His Lys Trp Glu Asn Phe Tyr Cys Glu Gln Gln Asp Pro Phe Val Cys Glu Ala Asp Cys Pro Ser Asp Trp Ser Ser Tyr Glu Gly His Cys Tyr Lys Pro Phe Asn Glu Pro Lys Asn Trp Ala Asp Ala Glu Asn Phe Cys Thr Lys Gln His Thr Gly Gly His Leu Val Ser Phe Gln Ser Thr Glu Glu Ala Asp Phe Val Val Lys Leu Ala Phe Gln Thr Phe Asp Tyr Gly Leu Phe Trp Phe Gly Leu Ser Lys Leu Trp Asn Gln Cys Asn Trp Gln Trp Ser Asn Ala Ala Met Leu Lys Tyr Thr Asp Trp Ala Glu Glu Ser Tyr Cys Val Tyr Phe Lys Ser Thr Asn Asn Lys Trp Arg Ser Leu Thr Cys Arg Met Leu Ala Asn Phe Val Cys Glu Phe Gln Ala'
    #
    # sequence = sequence.split(' ')
    #
    # aa = [key[i] for i in sequence]
    # print (''.join(aa))
    #
    for i, v in smiles.items():

        print ("'%s': '%s'," % (i, Chem.MolToSmarts(Chem.MolFromSmiles(v))))

    for i, v in smiles.items():

        molecule = Chem.MolFromSmiles(v)
        bit_string = AllChem.GetMorganFingerprintAsBitVect(
            molecule,
            2,
            nBits=512
        ).ToBitString()

        print ("'%s': '%s'," % (i, bit_string))
