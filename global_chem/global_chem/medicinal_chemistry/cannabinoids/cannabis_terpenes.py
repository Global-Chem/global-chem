#!/usr/bin/env python3
#
# GlobalChem - CannabisTerpenes
#
# -----------------------------

class CannabisTerpenes(object):

  def __init__(self):

    self.name = 'cannabis_terpenes'

  @staticmethod
  def get_smiles():

    smiles = {
      'borneol': 'CC1(C)C2CCC1(C)C(O)C2',
      'bornyl acetate': 'CC(=O)OC1CC2CCC1(C)C2(C)C',
      'camphene': 'CC1(C)C2CCC(C2)C1=C',
      'camphenehydrate': 'CC1(C2CCC(C2)C1(C)O)C',
      'camphor': 'CC1(C)C2CCC1(C)C(=O)C2',
      'delta-3 carene': 'CC1=CCC2C(C1)C2(C)C',
      'delta-4 carene': 'CC1CC2C(C2(C)C)C=C1',
      'carvacrol': 'CC(C)c1ccc(C)c(O)c1',
      'carvone': 'CC(=C)C1CC=C(C)C(=O)C1',
      'beta-cyclocitral': 'CC1=C(C=O)C(C)(C)CCC1',
      '1,4-cineol': 'CC(C)C12CCC(C)(CC1)O2',
      '1,8-cineol': 'CC12CCC(CC1)C(C)(C)O2',
      'citral b': 'CC(C)=CCC\C(C)=C/C=O',
      'citronellol': 'CC(CCO)CCC=C(C)C',
      'para cymene': 'CC(C)c1ccc(C)cc1',
      'para cymene-8-ol': 'CC1=CC=C(C=C1)C(C)(C)O',
      'dihydrocarveyl acetate': 'CC1CCC(CC1OC(C)=O)C(C)=C',
      'dihydrocarvone': 'CC1CCC(CC1=O)C(C)=C',
      'fenchyl alcohol': 'CC1(C)C2CCC(C)(C2)C1O',
      'fenchone': 'CC1(C)C2CCC(C)(C2)C1=O',
      'geraniol': r'CC(C)=CCC\C(C)=C\CO',
      'geranyl acetone': r'CC(C)=CCCC(/C)=C/CCC(C)=O',
      'limonene': 'CC(=C)C1CCC(=CC1)C',
      'linalool': 'CC(C)=CCCC(C)(O)C=C',
      'linalool oxide': 'CC(C)=CCCC(C)(O)C1CO1',
      'meta mentha-1,8-dien-5-ol': '',
      '1-methyl-4-iso-propenylbenzene': 'CC(=C)c1ccc(C)cc1',
      'myrcene': 'CC(C)=CCCC(=C)C=C',
      'nerol':r'CC(C)=CCC\C(C)=C/CO',
      'nerolidol': 'CC(C)=CCCC(C)=CCCC(C)(O)C=C',
      'beta-ocimene': 'CC(C)=CCC=C(C)C=C',
      'perillene': 'CC(C)=CCCc1cocc1',
      'alpha-phellandrene': 'CC(C)C1CC=C(C)C=C1',
      'beta-phellandrene': 'CC(C)C1CCC(=C)C=C1',
      '3-phenyl-2-methyl-prop-1-ene': 'CC(=C)Cc1ccccc1',
      'alpha-pinene': 'CC1=CCC2CC1C2(C)C',
      'beta-pinene': 'CC1(C)C2CCC(=C)C1C2',
      'alpha-pinene oxide': 'CC1(C)C2CC3OC3(C)C1C2',
      'pinocarveol': 'CC1(C)C2CC(O)C(=C)C1C2',
      'pinocarvone': 'CC1(C)C2CC1C(=C)C(=O)C2',
      'piperitenone': 'CC(C)=C1CCC(=CC1=O)C',
      'piperitone oxide': 'CC(C)C1CCC2(C)OC2C1=O',
      'piperitenone oxide': 'CC(C)=C1CCC2(C)OC2C1=O',
      'pulegone': 'CC1CCC(=C(C)C)C(=O)C1',
      'sabinene': 'CC(C)C12CCC(=C)C1C2',
      'sabiene hydrate': 'CC(C)C12CCC(C1C2)(C)O',
      'sabinol': 'CC(C)C12CC(O)C(=C)C1C2',
      'safranal': 'CC1=C(C=O)C(C)(C)CC=C1',
      'alpha thujene': 'CC1=CCC2(C1C2)C(C)C',
      'alpha terpinene': 'CC1=CC=C(CC1)C(C)C',
      'gamma terpinene': 'CC1=CCC(=CC1)C(C)C',
      'alpha terpiene-4-ol': 'CC1=CCC(CC1)(C(C)C)O',
      'alpha terpinolene': 'CC1=CCC(=C(C)C)CC1',
      'alpha terpineol': 'CC1=CCC(CC1)C(C)(C)O',
      'beta terpineol': 'CC(=C)C1CCC(CC1)(C)O',
      'thujyl alcohol': 'CC(C)C12CC1C(C)CC2O',
      'allo-aromadendrene': 'CC1CCC2C1C3C(C3(C)C)CCC2=C',
      'alpha-bergamotene': 'CC(C)=CCCC1(C)C2CC=C(C)C1C2',
      'beta bisabolene': 'CC1=CCC(CC1)C(=C)CCC=C(C)C',
      'alpha bisabolol': 'CC1=CCC(CC1)C(C)(CCC=C(C)C)O',
      'calamenene': 'CC(C)C1CCC(C)c2ccc(C)cc12',
      'caryophyllene': r'C\C1=C/CCC(=C)C2CC(C)(C)C2CC1',
      'alpha caryophyllene': 'CC1=CCC(C=CCC(=CCC1)C)(C)C',
      'beta caryophyllene': 'CC1=CCCC(=C)C2CC(C2CC1)(C)C',
      'alpha caryophyllene alcohol': 'CC1(CC2C(C1)C3(CCCC2(C3O)C)C)C',
      'isocarophyllene': 'CC1=CCCC(=C)C2CC(C2CC1)(C)C',
      'caryophyllene oxide': 'CC1(C)CC2C1CCC3(C)OC3CCC2=C',
      'alpha cedrene': 'CC1CCC2C13CC=C(C(C3)C2(C)C)C',
      'gamma cadinene': 'CC1=CC2C(CC1)C(=C)CCC2C(C)C',
      'delta cadinene': 'CC1=CC2C(CCC(=C2CC1)C)C(C)C',
      'alpha copaene': 'CC1=CCC2C3C1C2(CCC3C(C)C)C',
      'alpha cubebene': 'CC1CCC(C2C13C2C(=CC3)C)C(C)C',
      'alpha curcumene': 'CC1=CC=C(C=C1)C(C)CCC=C(C)C',
      'beta cucumene': 'CC1=CCC(=CC1)C(C)CCC=C(C)C',
      'gamma elemene': 'CC(=C1CCC(C(C1)C(=C)C)(C)C=C)C',
      'gamma eudesmol': 'CC1=C2CC(CCC2(CCC1)C)C(C)(C)O',
      'beta farnesene': 'CC(=CCCC(=CCCC(=C)C=C)C)C',
      'z beta farnesene': 'CCC(=C)CCC=C(C)CCC=C(C)C',
      'alpha farnesene': 'CC(=CCCC(=CCC=C(C)C=C)C)C',
      'farnesol': 'CC(C)=CCCC(C)=CCCC(C)=CCO',
      'farnesyl acetone': 'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=O',
      'alpha gurjunene': 'CC1CCC2C(C2(C)C)C3=C(CCC13)C',
      'guaiol': 'CC1CCC(CC2=C1CCC2C)C(C)(C)O',
      'beta humulene': 'CC1=CCC(C=CCC(=C)CCC1)(C)C',
      'humelene epoxide I': 'CC1=CCCC2(C(O2)CC(C=CC1)(C)C)C',
      'humelene epoxide II': 'CC1=CCC(C=CCC2(C(O2)CC1)C)(C)C',
      'ledol': 'CC1CCC2C1C3C(CCC2(C)O)C3(C)C',
      'longifolene': 'CC1(C)CCCC2(C)C3CCC(C13)C2=C',
      'epi-beta-sanalene': 'CC(=CCCC1(C2CCC(C2)C1=C)C)C',
      'alpha selinene': 'CC1=CCCC2(C1CC(CC2)C(=C)C)C',
      'beta selinene': 'CC(=C)C1CCC2(CCCC(=C)C2C1)C',
      'selina-3,7-diene': 'CC(C)C1=CC[C@@]2(C)CCC=C(C)C2C1',
      'selina-4,7-diene': 'CC1=C2CC(=CCC2(CCC1)C)C(C)C',
      'friedelin': 'CC1C(=O)CCC2C1(C)CCC3C2(C)CCC4(C)C5CC(C)(C)CCC5(C)CCC34C',
      'epifriedelanol': 'CC1C(O)CCC2C1(C)CCC3C2(C)CCC4(C)C5CC(C)(C)CCC5(C)CCC34C',
      'vomifoliol': r'CC(O)/C=C/C1(O)C(=CC(=O)CC1(C)C)C',
      'dihydrovomifoliol': 'CC1(C)CCCC2(C)OC(=O)C=C12',
      'beta ionone': 'CC1=C(C(CCC1)(C)C)C=CC(=O)C',
      'dihydroactinidiolide': 'CC1(CCCC2(C1=CC(=O)O2)C)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
