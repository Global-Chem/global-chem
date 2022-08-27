
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
    'ziconotide': 'CC1C(=O)NC(C(=O)NC2CSSCC3C(=O)NC(C(=O)NC(C(=O)NCC(=O)NC(C(=O)NC(CSSCC(C(=O)NC(CSSCC(C(=O)NC(C(=O)NCC(=O)NC(C(=O)NCC(=O)N1)CCCCN)CCCCN)N)C(=O)NC(C(=O)NCC(=O)NC(C(=O)N3)CO)C(C)O)NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC2=O)CO)CCCNC(=N)N)CC(C)C)CCSC)CC4=CC=C(C=C4)O)CC(=O)O)C(=O)N)CCCCN)CO)CCCNC(=N)N)CCCCN',
    'captopril': 'CC(CS)C(=O)N1CCCC1C(=O)O',
    'tirofiban': 'CCCCS(=O)(=O)NC(CC1=CC=C(C=C1)OCCCCC2CCNCC2)C(=O)O',
    'eptifibatide': 'C1CC2C(=O)NC(CSSCCC(=O)NC(C(=O)NCC(=O)NC(C(=O)NC(C(=O)N2C1)CC3=CNC4=CC=CC=C43)CC(=O)O)CCCCN=C(N)N)C(=O)N',
    'batroxobin': 'NC(CCSC)C(NC(C(C)C)C(NC(CC(C)C)C(NC(C(CC)([H])C)C(NC(CCCCNC(N)=N)C(NC(C(C)C)C(NC(C(CC)([H])C)C(NC(C)C(NC(CC(N)=O)C(NC(CC(C)C)C(NC(CC(C)C)C(NC(C(CC)([H])C)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C(C)C)C(NC(CO)C(NC(CC1=CC=C(O)C=C1)C(NC(C)C(NC(CCC(N)=O)C(NC(CCCCN)C(NC(CO)C(NC(CO)C(NC(CCC(O)=O)C(NC(CC(C)C)C(NC(C(C)C)C(NC(C(CC)([H])C)C(NC([H])C(NC([H])C(NC(CC(O)=O)C(NC(CCC(O)=O)C(NC(CS)C(NC(CC(O)=O)C(NC(C(CC)([H])C)C(NC(CC(N)=O)C(NC(CCC(O)=O)C(NC(CC1=CNC=N1)C(NC(C2CCCN2)C(NC(CC1=CC=CC=C1)C(NC(CC(C)C)C(NC(C)C(NC(CC1=CC=CC=C1)C(NC(CCSC)C(NC(CC1=CC=C(O)C=C1)C(NC(CC1=CC=C(O)C=C1)C(NC(CO)C(NC(C2CCCN2)C(NC(CCCCNC(N)=N)C(NC(CC1=CC=C(O)C=C1)C(NC(CC1=CC=CC=C1)C(NC(CS)C(NC([H])C(NC(CCSC)C(NC(C(C)([H])O)C(NC(CC(C)C)C(NC(C(CC)([H])C)C(NC(CC(N)=O)C(NC(CCC(N)=O)C(NC(CCC(O)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(C(C)C)C(NC(CC(C)C)C(NC(C(C)([H])O)C(NC(C)C(NC(C)C(NC(CC1=CNC=N1)C(NC(CS)C(NC(CC(N)=O)C(NC(CCCCNC(N)=N)C(NC(CCCCNC(N)=N)C(NC(CC1=CC=CC=C1)C(NC(CCSC)C(NC(CCCCNC(N)=N)C(NC(C(CC)([H])C)C(NC(CC1=CNC=N1)C(NC(CC(C)C)C(NC([H])C(NC(CCCCN)C(NC(CC1=CNC=N1)C(NC(C)C(NC([H])C(NC(CO)C(NC(C(C)C)C(NC(C)C(NC(CC(N)=O)C(NC(CC1=CC=C(O)C=C1)C(NC(CC(O)=O)C(NC(CCC(O)=O)C(NC(C(C)C)C(NC(C(C)C)C(NC(CCCCNC(N)=N)C(NC(CC1=CC=C(O)C=C1)C(NC(C2CCCN2)C(NC(CCCCN)C(NC(CCC(O)=O)C(NC(CCCCN)C(NC(CC1=CC=CC=C1)C(NC(C(CC)([H])C)C(NC(CS)C(NC(C2CCCN2)C(NC(CC(N)=O)C(NC(CCCCN)C(NC(CCCCN)C(NC(CCCCN)C(NC(CC(N)=O)C(NC(C(C)C)C(NC(C(CC)([H])C)C(NC(C(C)([H])O)C(NC(CC(O)=O)C(NC(CCCCN)C(NC(CC(O)=O)C(NC(C(CC)([H])C)C(NC(CCSC)C(NC(CC(C)C)C(NC(C(CC)([H])C)C(NC(CCCCNC(N)=N)C(NC(CC(C)C)C(NC(CC(O)=O)C(NC(CCCCNC(N)=N)C(NC(C2CCCN2)C(NC(C(C)C)C(NC(CCCCN)C(NC(CC(N)=O)C(NC(CO)C(NC(CCC(O)=O)C(NC(CC1=CNC=N1)C(NC(C(CC)([H])C)C(NC(C)C(NC(C2CCCN2)C(NC(CC(C)C)C(NC(CO)C(NC(CC(C)C)C(NC(C2CCCN2)C(NC(CO)C(NC(CC(N)=O)C(NC(C2CCCN2)C(NC(C2CCCN2)C(NC(CO)C(NC(C(C)C)C(NC([H])C(NC(CO)C(NC(C(C)C)C(NC(CS)C(NC(CCCCNC(N)=N)C(NC(C(CC)([H])C)C(NC(CCSC)C(NC([H])C(NC(CCC1=CNC2=C1C=CC=C2)C(NC([H])C(NC(C)C(NC(C(CC)([H])C)C(NC(C(C)([H])O)C(NC(C(C)([H])O)C(NC(CO)C(NC(CCC(O)=O)C(NC(CC(O)=O)C(NC(C(C)([H])O)C(NC(CC1=CC=C(O)C=C1)C(NC(C2CCCN2)C(NC(CC(O)=O)C(NC(C(C)C)C(NC(C2CCCN2)C(NC(CC1=CNC=N1)C(NC(CS)C(NC(C)C(NC(CC(N)=O)C(NC(C(CC)([H])C)C(NC(CC(N)=O)C(NC(CC(C)C)C(NC(CC1=CC=CC=C1)C(NC(CC(N)=O)C(NC(CC(N)=O)C(NC(C(C)([H])O)C(NC(C(C)C)C(NC(CS)C(NC(CCCCNC(N)=N)C(NC(CCC(O)=O)C(NC(C)C(NC(CC1=CC=C(O)C=C1)C(NC(CC(N)=O)C(NC([H])C(NC(CC(C)C)C(NC(C2CCCN2)C(NC(C)C(NC(CCCCN)C(NC(C(C)([H])O)C(NC(CC(C)C)C(NC(CS)C(NC(C)C(NC([H])C(NC(C(C)C)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC([H])C(NC([H])C(NC(C(CC)([H])C)C(NC(CC(O)=O)C(NC(C(C)([H])O)C(NC(CS)C(NC([H])C(NC([H])C(NC(CC(O)=O)C(NC(CO)C(NC([H])C(NC([H])C(NC(C2CCCN2)C(NC(CC(C)C)C(NC(C(CC)([H])C)C(NC(CS)C(NC(CC(N)=O)C(NC([H])C(NC(CCC(N)=O)C(NC(CC1=CC=CC=C1)C(NC(CCC(N)=O)C(NC([H])C(NC(C(CC)([H])C)C(NC(CC(C)C)C(NC(CO)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC([H])C(NC(CO)C(NC(CC(O)=O)C(NC(C2CCCN2)C(NC(CS)C(NC(C)C(NC(CCC(O)=O)C(NC(C2CCCN2)C(NC(CCCCNC(N)=N)C(NC(CCCCN)C(NC(C2CCCN2)C(NC(C)C(NC(CC1=CC=CC=C1)C(NC(CC1=CC=C(O)C=C1)C(NC(C(C)([H])O)C(NC(CCCCN)C(NC(C(C)C)C(NC(CC1=CC=CC=C1)C(NC(CC(O)=O)C(NC(CC1=CC=C(O)C=C1)C(NC(CC(C)C)C(NC(C2CCCN2)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(C(CC)([H])C)C(NC(CCC(N)=O)C(NC(CO)C(NC(C(CC)([H])C)C(NC(C(CC)([H])C)C(NC(C)C(NC([H])C(NC(CC(N)=O)C(NC(CCCCN)C(NC(C(C)([H])O)C(NC(C)C(NC(C(C)([H])O)C(NC(CS)C(NC(C2CCCN2)C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O',
    'hemocoagulase': 'NC(CC(O)=O)C(NC(CS)C(NC(CO)C(NC(CO)C(NC([H])C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CO)C(NC(CO)C(NC(CC1=CC=C(O)C=C1)C(NC(CCC(O)=O)C(NC([H])C(NC(CC1=CNC=N1)C(NC(CS)C(NC(CC1=CC=C(O)C=C1)C(NC(CCCCN)C(NC(C(C)C)C(NC(CC1=CC=CC=C1)C(NC(CCCCN)C(NC(CCC(N)=O)C(NC(CO)C(NC(CCCCN)C(NC(C(C)([H])O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(C)C(NC(CC(O)=O)C(NC(C)C(NC(CCC(O)=O)C(NC(CO)C(NC(CC1=CC=CC=C1)C(NC(CS)C(NC(C(C)([H])O)C(NC(CCCCN)C(NC(CCC(N)=O)C(NC(C(C)C)C(NC(CC(N)=O)C(NC([H])C(NC([H])C(NC(CC1=CNC=N1)C(NC(CC(C)C)C(NC(C(C)C)C(NC(CO)C(NC(CC(C)C)C(NC(CCC(O)=O)C(NC(CO)C(NC(CO)C(NC([H])C(NC(CCC(O)=O)C(NC(C)C(NC(CC(O)=O)C(NC(CC1=CC=CC=C1)C(NC(C(C)C)C(NC([H])C(NC(CCC(N)=O)C(NC(CC(C)C)C(NC(CC(C)C)C(NC(C)C(NC(CCC(N)=O)C(NC(CCCCN)C(NC(CC(C)C)C(NC(CCCCN)C(NC(CO)C(NC(C)C(NC(CCCCN)C(NC(CC(C)C)C(NC(CC1=CNC=N1)C(NC(C(C)C)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CC(C)C)C(NC([H])C(NC(CC(C)C)C(NC(CCCCNC(N)=N)C(NC(C)C(NC(CCC(N)=O)C(NC(CC(N)=O)C(NC(CCCCN)C(NC(CCC(O)=O)C(NC(CCCCN)C(NC(CCC(N)=O)C(NC(CS)C(NC(CO)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CO)C(NC(CC(O)=O)C(NC([H])C(NC(CO)C(NC(CO)C(NC(CC(C)C)C(NC(CO)C(NC(CC1=CC=C(O)C=C1)C(NC(CCC(O)=O)C(NC(CC(N)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CC(C)C)C(NC(CCC(O)=O)C(NC(CCC(O)=O)C(NC(CCC(O)=O)C(NC(CO)C(NC(CCCCN)C(NC(CCCCN)C(NC(CS)C(NC(CC(C)C)C(NC([H])C(NC(C(C)C)C(NC(CC1=CNC=N1)C(NC(CC(C)C)C(NC(CCC(O)=O)C(NC(C(C)([H])O)C(NC([H])C(NC(CC1=CC=CC=C1)C(NC(CC1=CNC=N1)C(NC(CCCCN)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CCC(O)=O)C(NC(CC(N)=O)C(NC(CC1=CC=CC=C1)C(NC(CC1=CC=C(O)C=C1)C(NC(CS)C(NC(CCC(O)=O)C(NC(CCC(N)=O)C(NC(CCC(N)=O)C(NC(CC(O)=O)C(NC(C2CCCN2)C(NC(CC1=CC=CC=C1)C(NC(C(C)C)C(NC(CS)C(NC(CCC(O)=O)C(NC(C)C(NC(CC(O)=O)C(NC(CS)C(NC(C2CCCN2)C(NC(CO)C(NC(CC(O)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CO)C(NC(CO)C(NC(CC1=CC=C(O)C=C1)C(NC(CCC(O)=O)C(NC([H])C(NC(CC1=CNC=N1)C(NC(CS)C(NC(CC1=CC=C(O)C=C1)C(NC(CCCCN)C(NC(C2CCCN2)C(NC(CC1=CC=CC=C1)C(NC(CC(N)=O)C(NC(CCC(O)=O)C(NC(C2CCCN2)C(NC(CCCCN)C(NC(CC(N)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(C)C(NC(CC(O)=O)C(NC(C)C(NC(CCC(O)=O)C(NC(CC(N)=O)C(NC(CC1=CC=CC=C1)C(NC(CS)C(NC(C(C)([H])O)C(NC(CCCCN)C(NC(CCC(N)=O)C(NC(CC1=CNC=N1)C(NC(C(C)([H])O)C(NC([H])C(NC([H])C(NC(CC1=CNC=N1)C(NC(CC(C)C)C(NC(C(C)C)C(NC(CO)C(NC(CC1=CC=CC=C1)C(NC(CCC(N)=O)C(NC(CO)C(NC(C(C)([H])O)C(NC(CCC(O)=O)C(NC(CCC(O)=O)C(NC(C)C(NC(CC(O)=O)C(NC(CC1=CC=CC=C1)C(NC(C(C)C)C(NC(C(C)C)C(NC(CCCCN)C(NC(CC(C)C)C(NC(C)C(NC(CC1=CC=CC=C1)C(NC(CCC(N)=O)C(NC(C(C)([H])O)C(NC(CC1=CC=CC=C1)C(NC(CC(O)=O)C(NC(CC1=CC=C(O)C=C1)C(NC([H])C(NC(CC(C)C)C(NC(CC1=CC=CC=C1)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CC1=CC=CC=C1)C(NC([H])C(NC(CC(C)C)C(NC(CO)C(NC(CCCCN)C(NC(CC(C)C)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CC(N)=O)C(NC(CCC(N)=O)C(NC(CS)C(NC(CC(N)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CCC(N)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CO)C(NC(CC(N)=O)C(NC(C)C(NC(C)C(NC(CCSC)C(NC(CC(C)C)C(NC(CCCCN)C(NC(CC1=CC=C(O)C=C1)C(NC(C(C)([H])O)C(NC(CC(O)=O)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(C)C(NC(CCC(O)=O)C(NC(CCC(O)=O)C(NC(CO)C(NC(CC1=CC=C(O)C=C1)C(NC(CS)C(NC(C(C)C)C(NC(CC1=CC=C(O)C=C1)C(NC(CC1=CC=CC=C1)C(NC(CCCCN)C(NC(CO)C(NC(C(C)([H])O)C(NC(CC(N)=O)C(NC(CC(N)=O)C(NC(CCCCN)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CCCCNC(N)=N)C(NC(CO)C(NC(CC(C)C)C(NC(C(C)([H])O)C(NC(CS)C(NC(CCCCNC(N)=N)C(NC(CCSC)C(NC(CC(C)C)C(NC(C)C(NC(CC(N)=O)C(NC(CC1=CC=CC=C1)C(NC(C(C)C)C(NC(CS)C(NC(CCC(O)=O)C(NC(CC1=CC=CC=C1)C(NC(CCC(N)=O)C(NC(C)C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O',
    'alpha cobrotoxin': 'NC(CC(C)C)C(NC(CCC(O)=O)C(NC(CS)C(NC(CC1=CNC=N1)C(NC(CC(N)=O)C(NC(CCC(N)=O)C(NC(CCC(N)=O)C(NC(CO)C(NC(CO)C(NC(CCC(N)=O)C(NC(C(C)([H])O)C(NC(C2CCCN2)C(NC(C(C)([H])O)C(NC(C(C)([H])O)C(NC(C(C)([H])O)C(NC([H])C(NC(CS)C(NC(CO)C(NC([H])C(NC([H])C(NC(CCC(O)=O)C(NC(C(C)([H])O)C(NC(CC(N)=O)C(NC(CS)C(NC(CC1=CC=C(O)C=C1)C(NC(CCCCN)C(NC(CCCCN)C(NC(CCCCNC(N)=N)C(NC(CCC1=CNC2=C1C=CC=C2)C(NC(CCCCNC(N)=N)C(NC(CC(O)=O)C(NC(CC1=CNC=N1)C(NC(CCCCNC(N)=N)C(NC([H])C(NC(CC1=CC=C(O)C=C1)C(NC(CCCCNC(N)=N)C(NC(C(C)([H])O)C(NC(CCC(O)=O)C(NC(CCCCNC(N)=N)C(NC([H])C(NC(CS)C(NC([H])C(NC(CS)C(NC(C2CCCN2)C(NC(CO)C(NC(C(C)C)C(NC(CCCCN)C(NC(CC(N)=O)C(NC([H])C(NC(C(CC)([H])C)C(NC(CCC(O)=O)C(NC(C(CC)([H])C)C(NC(CC(N)=O)C(NC(CS)C(NC(CS)C(NC(C(C)([H])O)C(NC(C(C)([H])O)C(NC(CC(O)=O)C(NC(CCCCNC(N)=N)C(NC(CS)C(NC(CC(N)=O)C(NC(CC(N)=O)C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O',
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
