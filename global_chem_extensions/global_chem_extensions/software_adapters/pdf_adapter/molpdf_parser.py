#!/usr/bin/env python3
#
# GlobalChemExtensions - MolPDF
#
# ------------------------------------


from molpdf import MolPDF, MolPDFParser

class MolPDFAdpater(object):

    __version__ = '0.0.1'


    def __init__(self,
                 smiles = [],
                 labels = [],
                 file_name = 'molecules.pdf',
                 include_failed_smiles = True,
                 title = 'CHEMICAL LIST BY MOLPDF',
                 ):

        self.title = title
        self.labels = labels
        self.smiles_list = smiles
        self.file_name = file_name
        self.include_failed_smiles = include_failed_smiles

        self.document = None

    def initialize_document(self):

        '''

        Initialize the Document

        '''

        self.document = MolPDF(name=self.file_name)
        self.document.add_title(self.title)
        self.document.add_spacer()

    def generate_document(self):

        '''

        Generate the Document

        '''

        self.initialize_document()

        self.document.generate(
            smiles=self.smiles_list,
            include_failed_smiles=self.include_failed_smiles
        )

    def parse_document(self):

        '''

        Parse a MolPDF Document

        '''

        self.document = MolPDFParser(self.file_name)
        smiles_list = self.document.extract_smiles()

        return smiles_list