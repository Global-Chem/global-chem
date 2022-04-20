# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

# imports
# -------

from global_chem import GlobalChem

def test_initializa_class():

    '''

    Test the intialize of the class, with parameters if it extends to that

    '''

    gc = GlobalChem()

def test_print_network():

    '''

    Test the printing of the network

    '''

    gc = GlobalChem()
    gc.print_globalchem_network()

def test_node_attribute():

    '''

    Test the Node Attribute to the GlobalChem Class

    '''

    gc = GlobalChem()
    nodes_list = gc.check_available_nodes()

    print (nodes_list)
    assert len(nodes_list) > 0
    assert 'emerging_perfluoroalkyls' in nodes_list

def test_build_global_chem_network():

    '''

    Test the Building of the GlobalChem Network

    '''

    gc = GlobalChem()
    gc.build_global_chem_network(print_output=True)
    molecules = gc.get_node_smiles('emerging_perfluoroalkyls')

    assert len(molecules) > 0

def test_fetch_all_globalchem_data():

    '''

    Test the fetching of all data

    '''

    gc = GlobalChem()

    all_smiles = gc.get_all_smiles()
    all_smarts = gc.get_all_smarts()
    all_names = gc.get_all_names()

    assert len(all_smiles) > 0
    assert len(all_smarts) > 0
    assert len(all_names) > 0

def test_fetch_node_data():

    '''

    Fetch the Node Data and Test

    '''

    gc = GlobalChem()

    gc.build_global_chem_network(print_output=True)

    smiles = gc.get_node_smiles('emerging_perfluoroalkyls')
    smarts = gc.get_node_smarts('emerging_perfluoroalkyls')

    assert len(smiles) > 0
    assert len(smarts) > 0

def test_get_and_set_node_values():

    '''

    Test the get and set of the nodes

    '''

    gc = GlobalChem()
    gc.build_global_chem_network(print_output=True, debugger=False)

    gc.set_node_value('emerging_perfluoroalkyls', {'some_data': ['bunny']})
    node_value = gc.get_node_value('emerging_perfluoroalkyls')

    assert len(node_value) > 0

def remove_node():

    '''

    Test the removal of a node

    '''

    gc = GlobalChem()
    gc.build_global_chem_network(print_output=True, debugger=False)

    gc.remove_node('emerging_perfluoroalkyls')
    keys = list(gc.network.keys())

    assert 'emerging_perfluoroalkyls' not in keys

def test_get_node():

    '''

    Test the get of a node

    '''

    gc = GlobalChem()
    gc.build_global_chem_network()
    node = gc.get_node('emerging_perfluoroalkyls')

    assert len(node) > 0

def test_fetch_smiles_by_iupac():

    '''

    Test the fetcher function for IUPAC by SMILES

    '''

    gc = GlobalChem()
    definition = gc.get_smiles_by_iupac(
        'benzene',
        return_network_path=False,          # Return the last found network path
        return_all_network_paths=True,      # Return all the found network paths
    )

    keys = list(definition.keys())

    assert 'C1=CC=CC=C1' in keys

def test_compute_common_score():

    '''

    Compute the Common Score of a molecule

    '''

    gc = GlobalChem()
    gc.build_global_chem_network(print_output=False, debugger=False)
    gc.compute_common_score('benzene', verbose=True)

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

    assert len(values) > 0

    # Test Fetching SMARTS

    values = gc.get_node_smarts('electrophilic_warheads_for_kinases')

    assert len(values) > 0

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
            'common_warheads_covalent_inhibitors',
            'privileged_scaffolds',
            'iupac_blue_book'
        ]
    )

    gc.print_deep_network()


