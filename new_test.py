class WaterBasedCoatings(object):

  @staticmethod
  def get_smiles():

      smiles = {
        'n,n-dimethylacrylamide': 'CN(C)C(=O)C=C',
        'acrylamide': 'C=CC(N)=O',
        'acrylic acid,': 'C=CC(=O)O',
        'methacrylic acid': 'CC(=C)C(=O)O',
        'styrene': 'C=Cc1ccccc1',
        'acrylonitrile': 'C=CC#N',
        'methyl methacrylate': 'C=CO',
        'butyl acrylate': 'CCCCOC(=O)C=C',
        'ethyl acrylate': 'CCOC(=O)C=C',
        '2-ethyl hexyl acrylate': 'CCCCC(CC)COC(=O)C=C',
        'hexamethylene diisocyanate': 'C(CCCN=C=O)CCN=C=O',            # Think there is a mistake here
        'toluene diisocyanate': 'CC1=CC=CC=C1.C(=[N-])=O.C(=[N-])=O',  # Think there is a mistake here
        'polyethylene glycol': 'COCCCCO',
        'biuret': 'C(=O)(N)NC(=O)N',
        'vinyl acetate': 'CC(=O)OC=C',
        'veova-10': 'CC(C)(C)CCCCCC(=O)OC=C',
      }

      return smiles

class SolventBasedCoatings(object):

  @staticmethod
  def get_smiles():
      smiles = {
        'phathlic anhydride': 'C1=CC=C2C(=C1)C(=O)OC2=O',
        'isopthalic acid': 'C1=CC(=CC(=C1)C(=O)O)C(=O)O',
        'terepthalic acid': 'C1=CC(=CC=C1C(=O)O)C(=O)O',
        'diethylene glycol': 'C(COCCO)O',
        'adipic acid': 'C(CCC(=O)O)CC(=O)O',
        'pentanedioyl dichloride': 'C(CC(=O)Cl)CC(=O)Cl',
        'benzene-1,4-diol': 'C1=CC(=CC=C1O)O',
        'ricinoleic acid': 'CCCCCCC(CC=CCCCCCCCC(=O)O)O',
        'oleic acid': r'CCCCCCCC/C=C\CCCCCCCC(=O)O',
        'linoleic acid': 'CCCCCC=CCC=CCCCCCCCC(=O)O',
        'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O',
        'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
      }

      return smiles

class ThermoPlastics(object):

  @staticmethod
  def get_smiles():
      smiles = {
        'bisphenol a': 'CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O',
        'bisphenol f': 'C1=CC(=CC=C1CC2=CC=C(C=C2)O)O',
        'chloromethyloxirane': 'C1C(O1)CCl',
        'formaldehyde': 'C=O',
        'epsilon caprolactone': 'C1CCC(=O)OCC1',
        'propylene oxide': 'CC1CO1',
        'ethylene oxide': 'C1CO1',
        'urea': 'C(=O)(N)N',
        'phenol': 'C1=CC=C(C=C1)O',
        'melamine': 'C1(=NC(=NC(=N1)N)N)N',
      }

      return smiles

class Rubbers(object):

    @staticmethod
    def get_smiles():
        smiles = {
          '2-chlorobuta-1,3-diene': 'C=CC(=C)Cl',
          'isoprene': 'CC(=C)C=C',
          'ethylidene norbornene': 'CC=C1CC2CC1C=C2',
          'dicyclopentadiene': 'C1C=CC2C1C3CC2C=C3',
          'vinyl norbornene': 'C=CC1CC2CC1C=C2',
          'vinylidene fluoride': 'C=C(F)F',
          'ethylene': 'C=C',
          'butadiene': 'C=CC=C',
        }

        return smiles

class NaturalFibers(object):

  @staticmethod
  def get_smiles():
      smiles = {
        'p-coumaric acid': 'C1=CC(=CC=C1C=CC(=O)O)O',
        'ferulic acid': 'COC1=C(C=CC(=C1)C=CC(=O)O)O',
        'syringyl alcohol': 'CC(=O)OCC1=CC(=C(C(=C1)OC)OC(=O)C)OC',
        'gallic acid': 'C1=C(C=C(C(=C1O)O)O)C(=O)O',
        'cellulose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
        'beta glucose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'starch': '', # What type of starch
        'alpha glucose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'hemicellulose': 'CC(C(C1=CC=CC=C1)O)NC.CC(C(C1=CC=CC=C1)O)NC',
        'arabinose': 'C1C(C(C(C(O1)O)O)O)O',
        'mannose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'glucose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'xylose': 'C1C(C(C(C(O1)O)O)O)O',
      }

      return smiles

class FluoroPolymers(object):

    @staticmethod
    def get_smiles():
        smiles = {
          'perfluorocycloalkene': 'C1(=C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)F)C(F)(F)F',
          'vinyl fluoride': 'C=CF',
          'vinylidene': 'C=C(F)F',
          'tetrafluoroethylene': 'C(=C(F)F)(F)F',
          'chlorotrifluoroethylene': 'C(=C(F)Cl)(F)F',
          'hexafluoropropylene': 'C(=C(F)F)(C(F)(F)F)F',
        }

        return smiles

class Silicones(object);

      @staticmethod
      def get_smiles():
        smiles = {
          'octamethylcyclotetrasiloxane': 'C[Si]1(O[Si](O[Si](O[Si](O1)(C)C)(C)C)(C)C)',
        }

        return smiles



