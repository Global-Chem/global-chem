from global_chem_package.global_chem.global_chem import GlobalChem

json_object = {}
file = open('test.json', 'w')

group = 10


if __name__ == '__main__':

    from global_chem import GlobalChem

    gc = GlobalChem()
    gc.build_global_chem_network(print_output=False, debugger=False)

    nodes = gc.get_all_nodes()

    # for node in nodes:
    #
    #     node_name = node.__name__
    #
    #     if node_name == 'Node' or node_name == 'GlobalChem' or node_name == 'CommonRegexPatterns':
    #         continue
    #
    #     smiles_list = node.get_smiles()
    #     iupac_names = smiles_list.keys()
    #
    #     for name in iupac_names:
    #
    #         template = '{ "id": "%s", "group": %s },' % (name, group)
    #         print (template)
    #
    #     group += 1

    for node in nodes[1:]:

        node_name = node.__name__

        if node_name == 'Node' or node_name == 'GlobalChem' or node_name == 'CommonRegexPatterns':
            continue

        smiles_list = node.get_smiles()
        iupac_names = smiles_list.keys()

        for name in iupac_names:
                template = '{ "source": "%s", "target": "%s" , "value": %s},' % (name, node_name, 1)
                print (template)


