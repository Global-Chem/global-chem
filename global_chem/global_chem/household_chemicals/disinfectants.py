#!/usr/bin/env python3
#
# GlobalChem - Disinfectants
#
# -------------------------------

class Disinfectants(object):

  def __init__(self):

    self.name = 'disinfectants'

  @staticmethod
  def get_smiles():

    smiles = {
'triclosan' : 'C1=CC(=C(C=C1Cl)O)OC2=C(C=C(C=C2)Cl)Cl',
'benzalkonium chloride' : '[Cl-].CCCCCCCCCCCCN+(C)CC1=CC=CC=C1',
'isopropyl alcohol' : 'CC(C)O',
'chlorohexidine' : r 'C1=CC(=CC=C1N/C(=N/C(=NCCCCCCN=C(/N=C(/NC2=CC=C(C=C2)Cl)\N)N)N)/N)Cl',
'hydrogen peroxide' : 'OO',
'potassium chlorate' : '[O-]Cl(=O)=O.[K+]',
'potassium iodate' : '[O-]I(=O)=O.[K+]',
'potassium permanganate' : '[O-]Mn(=O)=O.[K+]',
'phenol' : 'C1=CC=C(C=C1)O',
'formaldehyde' : 'C=O',
'boric acid' : 'B(O)(O)O',
'triclosan' : 'C1=CC(=C(C=C1Cl)O)OC2=C(C=C(C=C2)Cl)Cl',
'oligo(2-(2-ethoxy)-ethoxyethyl)guanidinium chloride' : 'Cl.NC(N)=N.NCCOCCOCCN',
'polyhexamethylene guanidine' : 'Cl.NC(N)=N.NCCCCCCN',
'hexadecyltrimethylammonium bromide' : ' CCCCCCCCCCCCCCCCN+(C)C.[Br-]',
'methyltrioctylammonium chloride' : 'CCCCCCCCN+(CCCCCCCC)CCCCCCCC.[Cl-]',
'methylbenzethonium chloride' : ' CC1=C(C=CC(=C1)C(C)(C)CC(C)(C)C)OCCOCCN+(C)CC2=CC=CC=C2.[Cl-]',
'tributyltetradecylphosphonium chloride' : 'CCCCCCCCCCCCCCP+(CCCC)CCCC.[Cl-]',
'didecyl dimethyl ammonium chloride' : 'CCCCCCCCCCN+(C)CCCCCCCCCC.[Cl-]',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {

 }

    return smarts
