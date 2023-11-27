#!/usr/bin/env python3
#
# GlobalChem - Mango Phytocompounds
# ----------------------------------

class MangoPhytocompounds(object):

  def __init__(self):

    self.name = 'mango_phytocompounds'

    @staticmethod
    def get_smiles():

        smiles = {
          'myricetin': 'C1=C(C=C(C(=C1O)O)O)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O',
          '9-cis-violaxanthin': 'C/C(=C\C=C\C=C(/C)\C=C\C=C(\C)/C=C/[C@]12[C@](O1)(C[C@H](CC2(C)C)O)C)/C=C/C=C(\C)/C=C/[C@]34[C@](O3)(C[C@H](CC4(C)C)O)C',
          'limonene': 'CC1=CCC(CC1)C(=C)C',
          'alpha-terpinolene': r"C/C(C)=C1CCC(C)C=C\1",
          'd-carvone': 'CC1=CC[C@@H](CC1=O)C(=C)C',
          'alpha phellandrene': '',
          'alpha humulene': '',
          'gamma terpinene': '',
          'alpha pinene': '',
          'trans caryophyllene': '',
          'sabinene': '',
          'carene': '',
          'cis-caryophyllene': '',
          'αlpha humulene': '',
          'germacrene d': '',
          'aromadendrene': '',
          'beta ubebene': '',
          'alpha cubebene': '',
          'alpha bourbonene': '',
          'beta elemene': '',

        }
        return smiles


