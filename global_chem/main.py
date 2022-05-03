from global_chem import GlobalChem


import selfies as sf

if __name__ == '__main__':

    gc = GlobalChem()

    gc.build_global_chem_network()

    nodes = [
        # 'pihkal',
        # 'schedule_one',
        # 'schedule_two',
        # 'schedule_three',
        # 'schedule_four',
        # 'schedule_five',
        # 'vitamins',
        # 'open_smiles',
        # 'common_warheads_covalent_inhibitors',
        # 'emerging_perfluoroalkyls',
        'amino_acids',
        # 'braf_inhibitors',
        # 'privileged_kinase_inhibitors',
    ]

    names_list = []
    selfies_list = []

    for node in nodes:

        names_list = list(gc.get_node_smiles(node).keys())
        smiles_list = list(gc.get_node_smiles(node).values())

        for i, smiles in enumerate(smiles_list):
            try:
                selfies = sf.encoder(smiles)
                names_list.append(names_list[i])
                selfies_list.append(selfies)
            except sf.EncoderError:
                pass

    print (selfies_list)
    print (names_list)

    with open('SELFIES_TO_IUPAC.txt', 'w') as f:


        for i, selfie in enumerate(selfies_list):

            f.write(selfie + '\t' + names_list[i] + '\n')