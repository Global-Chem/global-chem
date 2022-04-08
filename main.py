
from rdkit import Chem

compounds = {
    'microcrystalline_cellulose': 'COC1OC(CO)C(OC2OC(CO)C(OC)C(O)C2O)C(O)C1O',
    'corn_starch': 'COC4C(CO)OC(OCC2OC(OC1C(CO)OC(C)C(O)C1O)C(O)C(O)C2OC3OC(CO)C(OC)C(O)C3O)C(O)C4O',
    'dibasic_calcium_phosphate': 'O.O.OP(=O)([O-])[O-].[Ca+2]',
    'lactose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
    'pregelatinized_starch': 'COC1C(O)C(O)C(OCC2OC(OC3C(O)C(O)C(C)OC3CO)C(O)C(O)C2OC2OC(CO)C(OC)C(O)C2O)OC1CO',
    'hydroxypropyl_methylcellulose': 'CC(COCC1C(C(C(C(O1)OC2C(OC(C(C2OCC(C)O)OCC(C)O)OCC(C)O)COCC(C)O)OCC(C)O)OCC(C)O)OCC(C)O)O.COCC1C(C(C(C(O1)OC2C(OC(C(C2OC)OC)OC)COC)OC)OC)OC',
    'sodium_starch_glycolate': 'NN1CC2N(C(C)=O)C(C1)CC2',
    'sodium_lauryl_sulfate': 'CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]',
    'povidone': 'C1CC(=O)N(C1)C(CP)P',
    'croscarmellose_sodium': 'CC(=O)O.C(C(C(C(C(C=O)O)O)O)O)O.[Na]',
    'colloidal_silicon_dioxide': 'O=[Si]=O',
    'crospovidone': 'C1CC(=O)NC1',
    'stearic_acid': 'CCCCCCCCCCCCCCCCCC(=O)O',
    'magnesium_stearate': 'CCCCCCCCCCCCCCCCCC(=O)[O-].CCCCCCCCCCCCCCCCCC(=O)[O-].[Mg+2]'
}

if __name__ == '__main__':

    for k, v in compounds.items():

        v = Chem.MolToSmarts(Chem.MolFromSmiles(v))

        print (f"'{k}': '{v}',")
