class MangoFlavonoids(object):

    def __init__(self):

        self.name = 'mango_flavonoids'

    @staticmethod
    def get_smiles():

        smiles = {
            'quercetin-3-O-galactoside': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O)O)O',
            'quercetin-3-O-glucoside': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O',
            'quercetin-3-O-xyloside': 'C1C(C(C(C(O1)OC2=C(OC3=CC(=CC(=C3C2=O)O)O)C4=CC(=C(C=C4)O)O)O)O)O',
            'magniferin': 'C1=C2C(=CC(=C1O)O)OC3=C(C2=O)C(=C(C(=C3)O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O',
            'cyanidin': 'C1=CC(=C(C=C1C2=[O+]C3=CC(=CC(=C3C=C2O)O)O)O)O',
            'delphinidin': 'C1=C(C=C(C(=C1O)O)O)C2=[O+]C3=CC(=CC(=C3C=C2O)O)O.[Cl-]',
            'pelargonidin': 'C1=CC(=CC=C1C2=[O+]C3=CC(=CC(=C3C=C2O)O)O)O',
            'catechin': 'C1[C@@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O',
            'apigenin': 'C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O ',
            'luteolin': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O',
            'kaempferol': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O',
        }

        return smiles
    
    @staticmethod
    def get_smarts():

      smarts = {
              
          }
      
      return smarts