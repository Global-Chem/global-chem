# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

# imports
# -------
from indigo import *
from rdkit import Chem
from global_chem import GlobalChem

def test_rdkit_passing():

    '''

    Test the global chem class through rdkit parser

    '''

    total_smiles = GlobalChem().get_all_smiles()

    passing_molecules = []

    for i in range(0, len(total_smiles)):

        try:
            mol = Chem.MolFromSmiles(total_smiles[i])
            if mol:
                passing_molecules.append(total_smiles[i])
        except Exception as e:
            pass

    print ('Total Molecules: %s ' % (len(total_smiles)))
    print ('Total RDKit Passing Molecules: %s ' % len(passing_molecules))
    print ('Failed Compounds: %s' % list(set(total_smiles) - set(passing_molecules)))

def test_indigo_passing():

    '''

    Test the GlobalChem Passing through the indigo parser.

    '''

    total_smiles = GlobalChem().get_all_smiles()

    indigo = Indigo()

    success_compounds = []
    failed_compounds = []

    for smiles in total_smiles:

        try:
            molecule = indigo.loadMolecule(smiles)
            success_compounds.append(smiles)
        except IndigoException as e:
            failed_compounds.append(smiles)

    print ('Total Molecules: %s ' % (len(total_smiles)))
    print ('Total Indigo Passing Molecules: %s ' % len(success_compounds))
    print ('Failed Compounds: %s' % len(failed_compounds))
    print ('Failed Compounds: %s' % failed_compounds)

def test_building_global_chem_network_and_data_access():


    '''

    Test the Global Chem Building Network.

    '''

    # Check Node Function

    gc = GlobalChem()
    nodes_list = gc.check_available_nodes()

    assert 'emerging_perfluoroalkyls' in nodes_list

    # See if the network builds

    gc.build_global_chem_network(print_output=True)

    # Check Node Function Capability of getting SMILES/SMARTS

    smiles = gc.get_node_smiles('emerging_perfluoroalkyls')
    smarts = gc.get_node_smarts('emerging_perfluoroalkyls')

    assert 'C(=O)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)O' in smiles
    assert '[#6](=[#8])(-[#6](-[#6](-[#6](-[#6](-[#6](-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])-[#8]' in smarts

    all_smiles = gc.get_all_smiles()
    all_smarts = gc.get_all_smarts()
    all_names = gc.get_all_names()

    # Test a series of SMILES/SMARTS/IUPAC

    assert 'C(=O)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)O' in all_smiles
    assert '[#6](=[#8])(-[#6](-[#6](-[#6](-[#6](-[#6](-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])(-[#9])-[#9])-[#8]' in all_smarts
    assert 'perfluorohexanoic acid' in all_names
    assert 'arginine' in all_names
    assert 'OCCC1=C(C)[N+](CC2=CN=C(C)N=C2N)=CS1'in smiles

    # Test Node Get and Set Value Functionality

    gc = GlobalChem()
    gc.build_global_chem_network(print_output=True, debugger=False)

    gc.set_node_value('emerging_perfluoroalkyls', {'some_data': ['bunny']})
    node_value = gc.get_node_value('emerging_perfluoroalkyls')

    assert 'some_data' in node_value

    # Test Remove Node Functionality

    gc.remove_node('emerging_perfluoroalkyls')

    node = gc.network.get('emerging_perfluoroalkyls', None)

    assert node == None

def test_building_independent_networks():

    '''

    Test the independent networks of user defined.

    '''

    # Test Initiating Network and Fetching SMILES

    gc = GlobalChem(verbose=False)
    gc.initiate_network()
    gc.add_node('global_chem', 'common_monomer_repeating_units')
    gc.add_node('common_monomer_repeating_units','electrophilic_warheads_for_kinases')
    values = gc.get_node_smiles('common_monomer_repeating_units')

    assert 'ClC1=CC=CC=C1C2=CC=C(C3=CC=CC=C3)C(Br)=C2' in values

    # Test Fetching SMARTS

    values = gc.get_node_smarts('electrophilic_warheads_for_kinases')

    assert '[#6H]-[#6]' in values

def test_deep_layer_networks():

    '''

    Test the deep layer network graphs

    '''


    gc = GlobalChem()
    gc.initiate_deep_layer_network()
    gc.add_deep_layer(
        [
            'emerging_perfluoroalkyls',
            'montmorillonite_adsorption',
            'common_monomer_repeating_units'
        ]
    )
    gc.add_deep_layer(
        [
            'common_warhead_covalent_inhibitors',
            'privileged_scaffolds',
            'iupac_blue_book'
        ]
    )

    gc.print_deep_network()

    assert 'emerging_perfluoroalkyls' in gc.network


