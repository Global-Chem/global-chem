from global_chem import GlobalChem

if __name__ == '__main__':


    gc = GlobalChem()
    definition = gc.get_smiles_by_iupac(
        '(4R,4aR,7S,7aR,12bS)-3-methyl-2,4,4a,7,7a,13-hexahydro-1H-4,12-methanobenzofuro[3,2-e]isoquinoline-7,9-diol',
        distance_tolerance=1,
        return_partial_definitions=False,
        reconstruct_smiles=True,
    )
    print (definition)