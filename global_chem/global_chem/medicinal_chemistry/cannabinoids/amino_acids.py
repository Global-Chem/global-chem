#!/usr/bin/env python3
#
# GlobalChem - Cannabis Amino Acids
#
# ---------------------------------

class CannabisAminoAcids(object):

  def __init__(self):

    self.name = 'cannabis_amino_acids'

  @staticmethod
  def get_smiles():

    smiles = {
      'alanine': 'CC(C(=O)O)N',
      'asparatic acid': 'C(C(C(=O)O)N)C(=O)O',
      'cystine': 'C(C(C(=O)O)N)SSCC(C(=O)O)N',
      'glutamic acid': 'CCC(O)=O',
      'glycine': 'C(C(=O)O)N',
      'serine': 'C(C(C(=O)O)N)O',
      'arginine': 'C(CC(C(=O)O)N)CN=C(N)N',
      'histidine': 'C1=C(NC=N1)CC(C(=O)O)N',
      'isoleucine': 'CCC(C)C(C(=O)O)N',
      'leucine': 'CC(C)CC(C(=O)O)N',
      'lysine': 'C(CCN)CC(C(=O)O)N',
      'methionine': 'CSCCC(C(=O)O)N',
      'phenylalanine': 'C1=CC=C(C=C1)CC(C(=O)O)N',
      'proline': 'C1CC(NC1)C(=O)O',
      'threonine': 'CC(C(C(=O)O)N)O',
      'tryptophan': 'C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N',
      'tyrosine': 'C1=CC(=CC=C1CC(C(=O)O)N)O',
      'valine': 'CC(C)C(C(=O)O)N',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
