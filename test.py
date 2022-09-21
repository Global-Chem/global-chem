
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
    'poe-20-sorbitan monooleate': 'CCCCCCCCC=CCCCCCCCC(=O)OCCOCC(C1C(CC(O1)OCCO)OCCO)OCCO',
    'poe-20-sorbitan monolaurate': 'CCCCCCCCCCCC(=O)OCCOCC(C1C(C(CO1)OCCO)OCCO)OCCO',
    'sorbitan monooleate': 'CCCCCCCCC=CCCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
    'sorbitan monolaurate': 'CCCCCCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
    'sorbitan monostearate': 'CCCCCCCCCCCCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
    'poloxamer 188': 'CCCOC(C)COCCO',
    'poloxamer 407': 'CC1CO1.C1CO1',
    'poe-35-castor oil': r'CCCCCCC(CC=CCCCCCCCOC(=O)OCCOCC(COCCOC(=O)OCCCCCCCC=CCC(CCCCCC)O)OCCOC(=O)OCCCCCCCC=CCC(CCCCCC)O)O',
    'peg-660-12-hydroxystearate': 'CCCCCCC(C)O.[CH2]CCCCCCCCCC(=O)OCCO',
    'tocopheryl-peg 1000-succinate': 'CC1=C(C(=C(C2=C1OC(CC2)(C)CCCC(C)CCCC(C)CCCC(C)C)C)OC(=O)CCC(=O)OCCO)C',
    'sucrose laurate': 'CCCCCCCCCCCC(=O)OC1(C(C(C(C(O1)CO)O)O)O)C2(C(C(C(O2)CO)O)O)CO',
    'sucrose palmitate': 'CCCCCCCCCCCCCCCC(=O)OC1(C(C(C(C(O1)CO)O)O)O)C2(C(C(C(O2)CO)O)O)CO',
    'caprylyl glucoside': r'CCCCCCCCCCOC1C(C(C(C(O1)CO)O)O)O',
    'polyetetramethyl butyl phenol ether': 'CC(C)(C)CC(C)(C)C1=CC=C(C=C1)OCCO',
    'soybean lecithin': 'CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCC=CCCCCC',
    'egg lecithin': 'CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCC=CCCCCC',
    'diolelyl phosphatidyl choline': 'CCCCCCCCC=CCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCCCCCCCC',
    'distearoyl phosphatidyl glycerol': 'CCCCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(CO)O)OC(=O)CCCCCCCCCCCCCCCCC',
    'dimyristoyl phosphotidyl choline': 'CCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC',
    'dioctyl sodium sulfosuccinate': 'CCCCCCCCOC(=O)CC(C(=O)OCCCCCCCC)S(=O)(=O)[O-].[Na+]',
    'ethanol': 'CCO',
    'benzyl alcohol': 'C1=CC=C(C=C1)CO',
    'isopropanol': 'CC(C)O',
    'phenethyl alcohol': 'C1=CC=C(C=C1)CCO',
    'propylene glycol': 'CC(CO)O',
    'glycerol': 'C(C(CO)O)O',
    'butylene glycol': 'CC(CCO)O',
    'peg 400': 'CC(C)CCCCCCCCCCCCCCC(=O)OCCO',
    'transcutol': 'CCOCCOCCO',
    'tetraglycol': 'C(COCCOCCOCCO)O',
    'n-methyl pyrrolidone': 'CN1CCCC1=O',
    '2-pyrrolidone': 'C1CC(=O)NC1',
    'sodium deoxycholate': 'CC(CCC(=O)[O-])C1CCC2C1(C(CC3C2CCC4C3(CCC(C4)O)C)O)C.[Na+]',
    'caprylic acid': 'CCCCCCCC(=O)O',
    'sodium caprylate': 'CCCCCCCC(=O)[O-].[Na+]',
    'potassium sorbate': 'CC=CC=CC(=O)[O-].[K+]'
}

if __name__ == '__main__':
    # sequence = 'Asp Cys Ser Ser Gly Trp Ser Ser Tyr Glu Gly His Cys Tyr Lys Val Phe Lys Gln Ser Lys Thr Trp Ala Asp Ala Glu Ser Phe Cys Thr Lys Gln Val Asn Gly Gly His Leu Val Ser Leu Glu Ser Ser Gly Glu Ala Asp Phe Val Gly Gln Leu Leu Ala Gln Lys Leu Lys Ser Ala Lys Leu His Val Trp Leu Gly Leu Arg Ala Gln Asn Lys Glu Lys Gln Cys Ser Leu Gln Trp Ser Asp Gly Ser Ser Leu Ser Tyr Glu Asn Trp Leu Glu Glu Glu Ser Lys Lys Cys Leu Gly Val His Leu Glu Thr Gly Phe His Lys Trp Glu Asn Phe Tyr Cys Glu Gln Gln Asp Pro Phe Val Cys Glu Ala Asp Cys Pro Ser Asp Trp Ser Ser Tyr Glu Gly His Cys Tyr Lys Pro Phe Asn Glu Pro Lys Asn Trp Ala Asp Ala Glu Asn Phe Cys Thr Lys Gln His Thr Gly Gly His Leu Val Ser Phe Gln Ser Thr Glu Glu Ala Asp Phe Val Val Lys Leu Ala Phe Gln Thr Phe Asp Tyr Gly Leu Phe Trp Phe Gly Leu Ser Lys Leu Trp Asn Gln Cys Asn Trp Gln Trp Ser Asn Ala Ala Met Leu Lys Tyr Thr Asp Trp Ala Glu Glu Ser Tyr Cys Val Tyr Phe Lys Ser Thr Asn Asn Lys Trp Arg Ser Leu Thr Cys Arg Met Leu Ala Asn Phe Val Cys Glu Phe Gln Ala'
    #
    # sequence = sequence.split(' ')
    #
    # aa = [key[i] for i in sequence]
    # print (''.join(aa))
    #
    # for i, v in smiles.items():
    #
    #     print ("'%s': '%s'," % (i, Chem.MolToSmarts(Chem.MolFromSmiles(v))))

    for i, v in smiles.items():

        molecule = Chem.MolFromSmiles(v)
        bit_string = AllChem.GetMorganFingerprintAsBitVect(
            molecule,
            2,
            nBits=512
        ).ToBitString()

        print ("'%s': '%s'," % (i, bit_string))



