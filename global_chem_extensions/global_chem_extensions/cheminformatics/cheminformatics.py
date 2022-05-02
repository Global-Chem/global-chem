#!/usr/bin/env python3
#
# GlobalChemExtensions - Cheminformatics Master Object
#
# ----------------------------------------------------

class ExtensionsError(Exception):

    __version_error_parser__ = "0.0.1"
    __allow_update__ = False

    '''
    
    Raise an Extension Error if something is wrong. 
    
    '''
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class ChemInformatics(object):

    def __init__(self):

        self.name = 'cheminformatics'

    @staticmethod
    def sunburst_chemical_list(smiles_list, save_file=False):

        '''


        Sunburst a chemical list

        Arguments:
            smiles_list (String): list of smiles strings to analyze
            save_file (Boolean): whether the user would like it as a file

        Returns:
            sunburster_object (Sunburster Object): Sunburster Object

        '''

        from global_chem_extensions.cheminformatics.applications.sunburster import Sunburster

        sunburster_object = Sunburster(smiles_list, save_file)

        return sunburster_object


    @staticmethod
    def node_pca_analysis(
            smiles_list,
            morgan_radius = 1,
            bit_representation = 512,
            number_of_clusters = 5,
            number_of_components = 0.95,
            random_state = 0,
            file_name = 'pca_analysis.html',
            save_file = False,
            return_mol_ids = False,
    ):

        '''

        Perform a pca analysis on a node within globalchem, can be extended to lists outside of the dedicated SMILES.

        Arguments:
            smiles_list (List): list of SMILES that the user wants to cluster
            morgan_radius (Int): Morgan Radius of the chemical environment
            bit_representation (Int): Length of the bit representation
            number_of_clusters (Int): Number of clusters the user would like to do
            number_of_components (Int): How many PCA vectors to analyze
            random_state (Int):
            file_name (String): file name the user would like to input
            save_file (Bool): Whether the user wants to display the plot or save it.
            return_mol_ids (Bool): Return the molecule IDS for the user to mine.

        '''

        from global_chem_extensions.cheminformatics.applications.node_pca_analysis import PCAAnalysis

        pca_analysis = PCAAnalysis(
            smiles_list,
            morgan_radius,
            bit_representation,
            number_of_clusters,
            number_of_components,
            random_state,
            file_name,
            save_file=save_file,
            return_mol_ids = return_mol_ids
        )

        return_mol_ids = pca_analysis.conduct_analysis()

        if return_mol_ids:
            return return_mol_ids

    @staticmethod
    def smiles_to_amino_acids(
            smiles_list
    ):

        '''

        Arguments:
            smiles_list (List): List of the SMILES

        Returns:
            converted_list (List): Converted list of the SMILES to the amino acid converters

        '''

        from global_chem_extensions.cheminformatics.applications.amino_acid_converter import AminoAcidConverter

        converter = AminoAcidConverter()
        converted_list = []

        for smiles in smiles_list:

            converted_list.append(
                converted_list.append(converter.convert_smiles_to_amino_acid_sequence(smiles))
            )

        return converted_list

    @staticmethod
    def amino_acids_to_smiles(
            sequence
    ):

        '''

        Arguments:
            amino_acid_list (List): List of the Amino Acids

        Returns:
            sequence (String): SMILES of the amino acids

        '''

        from global_chem_extensions.cheminformatics.applications.amino_acid_converter import AminoAcidConverter

        converter = AminoAcidConverter()

        sequence = converter.convert_amino_acid_sequence_to_smiles(sequence)

        return sequence

    @staticmethod
    def filter_smiles_by_criteria(
            smiles_list,
            lipinski_rule_of_5=False,
            ghose=False,
            veber=False,
            rule_of_3=False,
            reos=False,
            drug_like=False,
            pass_all_filters=False
    ):

        '''

        Arguments:
            lipinski_rule_of_5 (Bool): Lipinski Rule of 5 Criteria
            ghose (Bool): Ghose Filter,
            veber (Bool): Veber Filter,
            rule_of_3 (Bool): Rule of 3 Filter
            reos (Bool): Reos Filter
            drug_like (Bool): Drug Like Filter
            pass_all_filters (Bool): whether the user would like to pass all the filters

        Returns:

            the filtered data set
        '''

        from global_chem_extensions.cheminformatics.applications.drug_design_filters import DrugDesignFilters

        drug_design_filters = DrugDesignFilters(
            smiles_list,
            lipinski_rule_of_5=lipinski_rule_of_5,
            ghose=ghose,
            veber=veber,
            rule_of_3=rule_of_3,
            reos=reos,
            drug_like=drug_like,
            pass_all_filters=pass_all_filters
        )

        return drug_design_filters.filter()

    @staticmethod
    def convert_to_networkx(network):

        '''

        Arguments:
            network (Dict): Convert to a networkx object for interoperability

        Returns:
            converted_network (Graph Object): Networkx Graph object

        '''

        from global_chem_extensions.cheminformatics.applications.networkx_adapter import NetworkxAdapter

        network_adapter = NetworkxAdapter()
        converted_network = network_adapter.convert(network)

        return converted_network

    @staticmethod
    def scatter_deep_layer_network(deep_layer_network, height=800, width=1700, save_file=False):

        '''

        Arguments:
            deep_layer_network (Object): DGN coming from GlobalChem

        '''

        from global_chem_extensions.cheminformatics.applications.deep_layer_scatter import DeepLayerScatter

        deep_layer_scatter = DeepLayerScatter(deep_layer_network, save_file=save_file)
        deep_layer_scatter.scatter(height=height, width=width)

    @staticmethod
    def smarts_pattern_identifier(host='0.0.0.0', port=5000, debugger=True):

        '''

        Launch a Flask App for the SMARTS Pattern Identifier

        '''

        from global_chem_extensions.cheminformatics.applications.smarts_pattern_identifier import SmartsPatternIdentifier

        spi = SmartsPatternIdentifier()
        spi.launch_app(host=host, port=port, debug=debugger)

    @staticmethod
    def find_protonation_states(
            smiles_list,
            min_ph=6.4,
            max_ph=8.4,
            pka_precision=1.0,
            max_variants=128,
            label_states=False
    ):

        '''

        Find the protonation states of a SMILES string

        Returns:
            states (Dict): States of the SMILES input.

        '''

        from global_chem_extensions.cheminformatics.applications.dimorphite_dl import DimorphiteAdapter

        dimorphite_adapter = DimorphiteAdapter(
            smiles_list,
            min_ph,
            max_ph,
            pka_precision,
            max_variants,
            label_states
        )

        states = dimorphite_adapter.run()

        return states

    @staticmethod
    def verify_smiles(
            smiles_list,
            partial_smiles=True,
            rdkit=False,
            pysmiles=False,
            molvs=False,
            deepsmiles=False,
            partial=False,
            selfies=False,
            return_failures=False
    ):

        '''

        Arguments:
            smiles_list (List): List of smiles
            rdkit (Bool): Whether the RDKit flag is true
            partial_smiles (Bool): Whether to pass the partial smiles validation
            pysmiles (Bool): Whether to pass the validation through pysmiles
            molvs (Bool): Whether to pass the validation through MolVS
            deepsmiles (Bool): deepSMILES validation for machine learning
            partial (Bool): whether the user would like to have partial fragments
            selfies (Bool): whether the user would like to pass in the selfies for machine learning
            return_failures (Bool): whether the user would like to have failures returned.

        '''

        from global_chem_extensions.cheminformatics.applications.partial_smiles import PartialSmilesValidation

        psv = PartialSmilesValidation(
            partial=partial
        )

        successes, failures = psv.validate(
            smiles_list,
            rdkit=rdkit,
            partial_smiles=partial_smiles,
            pysmiles=pysmiles,
            molvs=molvs,
            deepsmiles=deepsmiles,
            selfies=selfies
        )

        if return_failures:
            return successes, failures
        else:
            return successes

    @staticmethod
    def smiles_to_pdf(
            smiles = [],
            labels = [],
            file_name = 'molecules.pdf',
            include_failed_smiles = True,
            title = 'CHEMICAL LIST BY MOLPDF',
    ):

        '''


        Arguments:
            smiles (List): List of smiles you want to parse
            labels (List): List of labels you want to parse
            file_name (String): File name you want the output to be
            include_failed_smiles (Bool): Whether to include SMILES that didn't render.
            title (String): Title of the PDF
        '''

        from global_chem_extensions.cheminformatics.applications.molpdf_parser import MolPDFAdapter

        molpdf_adapter = MolPDFAdapter(
            smiles = smiles,
            labels = labels,
            file_name = file_name,
            include_failed_smiles = include_failed_smiles,
            title = title
        )

        molpdf_adapter.generate_document()

    @staticmethod
    def pdf_to_smiles(
            file_name
    ):

        '''

        Arguments:
             file_name (String): File name for the pdf to be parsed

        Returns:
            molecules (List): List of molecules

        '''

        from global_chem_extensions.cheminformatics.applications.molpdf_parser import MolPDFAdapter

        molpdf_adapter = MolPDFAdapter(
            file_name = file_name
        )

        molecules = molpdf_adapter.parse_document()

        return molecules

    @staticmethod
    def encode_smiles(
            smiles_list,
            max_length = 120
    ):

        '''

        Arguments:

            smiles_list (List): List of SMILES
            max_length (Int): List of the encoded SMILES

        Returns:
            encoded_list (List): List of the encoded SMILES

        '''

        from global_chem_extensions.cheminformatics.applications.one_hot_encoding import SmilesOneHotEncoder

        encoder = SmilesOneHotEncoder(
            smiles_list = smiles_list,
            max_length = max_length
        )

        encoded_list = encoder.encode()

        return encoded_list

    @staticmethod
    def decode_smiles(
            smiles_list
    ):

        '''

        Arguments:

            smiles_list (List): List of SMILES

        Returns:
            decoded_list (List): List of the decoded SMILES

        '''

        from global_chem_extensions.cheminformatics.applications.one_hot_encoding import SmilesOneHotEncoder

        encoder = SmilesOneHotEncoder(
            smiles_list = smiles_list,
        )

        decoded_list = encoder.decode()

        return decoded_list

    @staticmethod
    def visualize_smarts(smarts):

        '''

        Arguments:
            smarts (String): Viusalize the SMARTS string

        '''

        from global_chem_extensions.cheminformatics.applications.smarts_visualizer import SmartsVisualizer

        visualizer = SmartsVisualizer(
            smarts
        )

        return visualizer.get_image()