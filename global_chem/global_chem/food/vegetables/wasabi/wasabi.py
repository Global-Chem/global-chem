#!/usr/bin/env python3
#
# GlobalChem - Wasabi
#
# -----------------------------------

class wasabi(object):

    def __init__(self):

        self.name = 'Wasabi'

    @staticmethod
    def get_smiles():


        smiles = {

        '6-O-caffeoylsucrose': r'C1=CC(=C(C=C1/C=C/C(=O)OCC2C(C(C(C(O2)OC3(C(C(C(O3)CO)O)O)CO)O)O)O)O)O',
        '6-O-feruloylsucrose': r'COC1=C(C=CC(=C1)/C=C/C(=O)OC[C@@H]2[C@H]([C@@H]([C@H]([C@H](O2)O[C@@]3([C@H]([C@@H]([C@H](O3)CO)O)O)CO)O)O)O)O',
        'ferulic Acid': r'COC1=C(C=CC(=C1)/C=C/C(=O)O)O',
        '6-methylsulfinyl-hexyl glucosinolate': 'CS(=O)CCCCCCC(=NOS(=O)(=O)O)S[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O',
        'isovitexin 4-O-glucoside': 'C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C=C(C(=C3O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O',
        'luteolin 7,3’-diglucoside': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O)O',
        'kaempferol 3-O-rutinoside': 'C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC=C(C=C5)O)O)O)O)O)O)O',
        'allyl isothiocyanate': 'C=CCN=C=S',
        '6-methylsulfinylhexyl isothiocyanate': 'CS(=O)CCCCCCN=C=S',
        'caffeic acid': r'C1=CC(=C(C=C1/C=C/C(=O)O)O)O',
        '1,2’-di-O-trans-sinapoyl gentiobiose' : r'O=C(O[C@@H]1[C@@H](O[H])[C@H](O[H])[C@@H](O[H])[C@H](CO[C@H]2[C@H](OC(/C=C/C3=CC(O[H])=C(O[H])C(O[H])=C3)=O)[C@@H](O[H])[C@H](O[H])[C@@H](CO[H])O2)O1)/C=C/C4=CC(O[H])=C(O[H])C(O[H])=C4',
        'apigenin 8-C-glucoside': 'C[Si](C (C)OCC1C(C(C(C(O1)C2=C(C=C(C3=C2OC(=CC3=O)C4=CC=C(C=C4)O[Si](C)(C)C)O[Si](C)(C)C)O[Si](C)(C)C)O[Si](C)(C)C)O[Si](C)(C)C)O[Si](C)(C)C',
        'p-coumaric acid': r'C1=CC(=CC=C1/C=C/C(=O)O)O',
        'n-butyl isothiocyanate': 'CCCCN=C=S',
        '3-butenyl isothiocyanate': 'C=CCCN=C=S',
        '4-pentenyl isothiocyanate': 'C=CCCCN=C=S',
        '5-hexenyl isothiocyanate': ' C=CCCCCN=C=S',
        'beta-Phenylethylisothiocyanate': 'C1=CC=C(C=C1)CCN=C=S',
        '7-methylthioheptyl isothiocyanate': 'CSCCCCCCCN=C=S', 
        '5-methylsulfinylpentyl isothiocyanate': 'CS(=O)CCCCCN=C=S',
        '7-methylsulfinylheptyl isothiocyanate': 'CS(=O)CCCCCCCN=C=S',
        'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
        'linolenic acid': r'CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O',
        'oleic acid': r'CCCCCCCC/C=C\CCCCCCCC(=O)O',
        'sinapinic acid': r'COC1=CC(=CC(=C1O)OC)/C=C/C(=O)O',
        '3,4-dimethoxy-trans-cinnamic acid': r'COC1=C(C=C(C=C1)/C=C/C(=O)O)OC',
        'sinigrin,2-propenyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CC=C)=N/OS(=O)(O)=O)C1O',
        '1-methylethyl glucisinolate': r'OC1C(O)C(CO)OC(S/C(C(C)C)=N/OS(=O)(O)=O)C1O',
        'glucpalyssin,5-methylsylfinylpentyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CCCCCS(C)=O)=N/OS(=O)(O)=O)C1O',
        'gluconapin,3-butenyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CCC=C)=N/OS(=O)(O)=O)C1O',
        'glucochlearin 1-methylpropyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(C(C)CC)=N/OS(=O)(O)=O)C1O',
        'glucohesperin, 6-methylsylfinylhexyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CCCCCCS(C)=O)=N/OS(=O)(O)=O)C1O',
        'glucoibarin 7-methylsylfinylheptyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CCCCCCCS(C)=O)=N/OS(=O)(O)=O)C1O',
        'glucobrassicanapin, 4-pentenyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CCCC=C)=N/OS(=O)(O)=O)C1O',
        '5-hexenyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CCCCC=C)=N/OS(=O)(O)=O)C1O',
        '4-methoxyglucobrassicin,4-methoxy-3-indolylmethyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CC2=CNC3=CC=CC(OC)=C23)=N/OS(=O)(O)=O)C1O',
        'neoglucobrassicin 1-methoxy-3indolylmethyl glucosinolate': r'OC1C(O)C(CO)OC(S/C(CC2=CN(OC)C3=CC=CC=C23)=N/OS(=O)(O)=O)C1O',
        '1-(3’’,4’’-dihydroxy-5’’-methoxy)-O-trans-cinnamoyl gentiobiose': r'O=C(/C=C/C1=CC(OC)=C(O)C(O)=C1)OC2OC(COC3C(O)C(O)C(O)C(CO)O3)C(O)C(O)C2O',
        '1,2’-di-(3’’,4’’-dihydroxy-5’’-methoxy)-O-trans-cinnamoyl gentiobiose': r'O=C(/C=C/C1=CC(OC)=C(O)C(O)=C1)OC2OC(COC3C(OC(/C=C/C4=CC(OC)=C(O)C(O)=C4)=O)C(O)C(O)C(CO)O3)C(O)C(O)C2O',
        '1-(3’’,4’’-dihydroxy-5’’-methoxy)-O-trans-cinnamoyl-2’-trans-sinapoyl gentiobiose': r'O=C(/C=C/C1=CC(OC)=C(O)C(O)=C1)OC2OC(COC3C(OC(/C=C/C4=CC(OC)=C(O)C(OC)=C4)=O)C(O)C(O)C(CO)O3)C(O)C(O)C2O',
        '1-O-trans-caffeoyl-2’-O-trans-sinapoyl gentiobiose': r'O=C(/C=C/C1=CC=C(O)C(O)=C1)OC2OC(COC3C(OC(/C=C/C4=CC(OC)=C(O)C(OC)=C4)=O)C(O)C(O)C(CO)O3)C(O)C(O)C2O',
        '1,2’-di-O-trans-sinapoyl gentiobiose': r'O=C(/C=C/C1=CC(OC)=C(O)C(OC)=C1)OC2OC(COC3C(OC(/C=C/C4=CC(OC)=C(O)C(OC)=C4)=O)C(O)C(O)C(CO)O3)C(O)C(O)C2O',
        '1-O-trans-feruoyl-2’O-trans-sinapoyl gentiobiose': r'O=C(/C=C/C1=CC(OC)=C(O)C(OC)=C1)OC2OC(COC3C(OC(/C=C/C4=CC(OC)=C(O)C=C4)=O)C(O)C(O)C(CO)O3)C(O)C(O)C2O',
        'isosaponarin': 'OCC1C(O)C(O)C(O)C(C(C(O)=C2)=C(O)C3=C2OC(C4=CC=C(OC5OC(O)C(O)C(O)C5CO)C=C4)=CC3=O)O1',
        'isovitexin': 'OCC1C(O)C(O)C(O)C(C(C(O)=C2)=C(O)C3=C2OC(C4=CC=C(O)C=C4)=CC3=O)O1',
        'isoorientin': 'OCC1C(O)C(O)C(O)C(C(C(O)=C2)=C(O)C3=C2OC(C4=CC=C(O)C(O)=C4)=CC3=O)O1',
        '7-O-trans-sinapoylisovitexin-4’-O-β-D-glucopyranoside': r'OCC1C(O)C(O)C(O)C(C(C(OC(/C=C/C2=CC(OC)=C(O)C(OC)=C2)=C)=C3)=C(O)C4=C3OC(C5=CC=C(OC6OC(O)C(O)C(O)C6CO)C=C5)=CC4=O)O1',

        
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
           
        }

        return smarts
