def _get_common_amino_acid_protecting_groups():

    alpha_amino_removed_by_acid_smiles = {
        'tert-butyloxycarbonyl': 'O=COC(C)(C)C',
        'trityl': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        '3,5-dimethoxyphenylisoproxycarbonyl': 'COC1=CC(C(C)(OC=O)C)=CC(OC)=C1',
        '2-(4-biphenyl)isopropoxycarbonyl': 'CC(C)(OC=O)C(C=C1)=CC=C1C2=CC=CC=C2',
        '2-nitrophenylsulfenyl': 'SC1=CC=CC=C1[N+]([O-])=O'
    }

    alpha_amino_removed_by_acid_smarts = {
        'tert-butyloxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        'trityl': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '3,5-dimethoxyphenylisoproxycarbonyl': '[#6]-[#8]-[#6]1:[#6]:[#6](-[#6](-[#6])(-[#8]-[#6]=[#8])-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:1',
        '2-(4-biphenyl)isopropoxycarbonyl': '[#6]-[#6](-[#6])(-[#8]-[#6]=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '2-nitrophenylsulfenyl': '[#16]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
    }

    alpha_amino_removed_by_acid_acronym_smiles = {
        'boc': 'O=COC(C)(C)C',
        'trt': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        'ddz': 'COC1=CC(C(C)(OC=O)C)=CC(OC)=C1',
        'bpoc': 'CC(C)(OC=O)C(C=C1)=CC=C1C2=CC=CC=C2',
        'nps': 'SC1=CC=CC=C1[N+]([O-])=O'
    }

    alpha_amino_removed_by_acid_acronym_smarts = {
        'boc': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        'trt': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'ddz': '[#6]-[#8]-[#6]1:[#6]:[#6](-[#6](-[#6])(-[#8]-[#6]=[#8])-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:1',
        'bpoc': '[#6]-[#6](-[#6])(-[#8]-[#6]=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'nps': '[#16]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
    }

    alpha_amino_removed_by_base_smiles = {
        '9-fluorenylmethoxycarbonyl': 'O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        '2-(4-nitrophenylsulfonyl)ethoxycarbonyl': 'O=COCCS(=O)(C1=CC=C([N+]([O-])=O)C=C1)=O',
        '(1,1-dioxobenzo[b]thiophene-2-yl)methyloxycarbonyl': 'O=COCC1=CC2=CC=CC=C2S1(=O)=O',
        '(1,1-dioxonaptho[1,2-b]thiophene-2-yl)methyloxycarbonyl': 'O=COCC1=CC2=CC=C3C=CC=CC3=C2S1(=O)=O',
        '1-(4,4-dimethyl-2,6-dioxocyclohex-1-ylidene)-3-methylbutyl': r'CC(C)C/C=C1C(CC(C)(C)CC\1=O)=O',
        '2,7-di-tert-butyl-fmoc': 'CC1=CC(C2COC=O)=C(C=C1)C3=C2C=C(C(C)(C)C)C=C3',
        '2-fluoro-fmoc': 'FC1=CC2=C(C(C=CC=C3)=C3C2COC=O)C=C1',
        '2-monoisooctyl-fmoc': 'O=COCC1C2=C(C=CC=C2)C3=C1C=C(C(CCCCCCC)=O)C=C3',
        '2,7-diisooctyl-fmoc': 'O=COCC1C2=C(C=CC(C(CCCCCCC)=O)=C2)C3=C1C=C(C(CCCCCCC)=O)C=C3',
        'tetrachlorophthaloyl': 'O=CC1=C(Cl)C(Cl)=C(Cl)C(Cl)=C1C=O',
        '2-[phenyl(methyl)sulfonio])ethyloxycarbonyltetrafluoroborate': 'C[S+](CCOC=O)C1=CC=CC=C1',
        'ethanesulfonylethoxycarbonyl': 'O=COC(S(=O)(CC)=O)C',
        '2-(4-sulfophenylsulfonyl)ethoxycarbonyl': 'O=COCCS(=O)(C1=CC=C(S(=O)(O)=O)C=C1)=O',
    }

    alpha_amino_removed_by_base_smarts = {
        '9-fluorenylmethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        '2-(4-nitrophenylsulfonyl)ethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]-[#16](=[#8])(-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1)=[#8]',
        '(1,1-dioxobenzo[b]thiophene-2-yl)methyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1=[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#16]-1(=[#8])=[#8]',
        '(1,1-dioxonaptho[1,2-b]thiophene-2-yl)methyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1=[#6]-[#6]2:[#6]:[#6]:[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3:[#6]:2-[#16]-1(=[#8])=[#8]',
        '1-(4,4-dimethyl-2,6-dioxocyclohex-1-ylidene)-3-methylbutyl': '[#6]-[#6](-[#6])-[#6]-[#6]=[#6]1-[#6](-[#6]-[#6](-[#6])(-[#6])-[#6]-[#6]-1=[#8])=[#8]',
        '2,7-di-tert-butyl-fmoc': '[#6]-[#6]1:[#6]:[#6]2-[#6](-[#6]-[#8]-[#6]=[#8])-[#6]3:[#6](-[#6]:2:[#6]:[#6]:1):[#6]:[#6]:[#6](:[#6]:3)-[#6](-[#6])(-[#6])-[#6]',
        '2-fluoro-fmoc': '[#9]-[#6]1:[#6]:[#6]2:[#6](-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3-[#6]-2-[#6]-[#8]-[#6]=[#8]):[#6]:[#6]:1',
        '2-monoisooctyl-fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]-1:[#6]:[#6](-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])=[#8]):[#6]:[#6]:2',
        '2,7-diisooctyl-fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](:[#6]:[#6]:[#6](-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])=[#8]):[#6]:2)-[#6]2:[#6]-1:[#6]:[#6](-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])=[#8]):[#6]:[#6]:2',
        'tetrachlorophthaloyl': '[#8]=[#6]-[#6]1:[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6]:1-[#6]=[#8]',
        '2-[phenyl(methyl)sulfonio])ethyloxycarbonyltetrafluoroborate': '[#6]-[#16+](-[#6]-[#6]-[#8]-[#6]=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'ethanesulfonylethoxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#16](=[#8])(-[#6]-[#6])=[#8])-[#6]',
        '2-(4-sulfophenylsulfonyl)ethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]-[#16](=[#8])(-[#6]1:[#6]:[#6]:[#6](-[#16](=[#8])(-[#8])=[#8]):[#6]:[#6]:1)=[#8]',
    }

    alpha_amino_removed_by_base_acronym_smiles = {
        'fmoc': 'O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        'nsc': 'O=COCCS(=O)(C1=CC=C([N+]([O-])=O)C=C1)=O',
        'bsmoc': 'O=COCC1=CC2=CC=CC=C2S1(=O)=O',
        'alpha-nsmoc': 'O=COCC1=CC2=CC=C3C=CC=CC3=C2S1(=O)=O',
        'ivdde':  r'CC(C)C/C=C1C(CC(C)(C)CC\1=O)=O',
        'fmoc*': 'CC1=CC(C2COC=O)=C(C=C1)C3=C2C=C(C(C)(C)C)C=C3',
        'fmoc(fmoc(2f))': 'FC1=CC2=C(C(C=CC=C3)=C3C2COC=O)C=C1',
        'mio-fmoc': 'O=COCC1C2=C(C=CC=C2)C3=C1C=C(C(CCCCCCC)=O)C=C3',
        'dio-fmoc': 'O=COCC1C2=C(C=CC(C(CCCCCCC)=O)=C2)C3=C1C=C(C(CCCCCCC)=O)C=C3',
        'tcp': 'O=CC1=C(Cl)C(Cl)=C(Cl)C(Cl)=C1C=O',
        'pms': 'C[S+](CCOC=O)C1=CC=CC=C1',
        'esc': 'O=COC(S(=O)(CC)=O)C',
        'sps': 'O=COCCS(=O)(C1=CC=C(S(=O)(O)=O)C=C1)=O',
    }

    alpha_amino_removed_by_base_acronym_smarts = {
        'fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        'nsc': '[#8]=[#6]-[#8]-[#6]-[#6]-[#16](=[#8])(-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1)=[#8]',
        'bsmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1=[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#16]-1(=[#8])=[#8]',
        'alpha-nsmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1=[#6]-[#6]2:[#6]:[#6]:[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3:[#6]:2-[#16]-1(=[#8])=[#8]',
        'ivdde': '[#6]-[#6](-[#6])-[#6]-[#6]=[#6]1-[#6](-[#6]-[#6](-[#6])(-[#6])-[#6]-[#6]-1=[#8])=[#8]',
        'fmoc*': '[#6]-[#6]1:[#6]:[#6]2-[#6](-[#6]-[#8]-[#6]=[#8])-[#6]3:[#6](-[#6]:2:[#6]:[#6]:1):[#6]:[#6]:[#6](:[#6]:3)-[#6](-[#6])(-[#6])-[#6]',
        'fmoc(fmoc(2f))': '[#9]-[#6]1:[#6]:[#6]2:[#6](-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3-[#6]-2-[#6]-[#8]-[#6]=[#8]):[#6]:[#6]:1',
        'mio-fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]-1:[#6]:[#6](-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])=[#8]):[#6]:[#6]:2',
        'dio-fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](:[#6]:[#6]:[#6](-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])=[#8]):[#6]:2)-[#6]2:[#6]-1:[#6]:[#6](-[#6](-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])=[#8]):[#6]:[#6]:2',
        'tcp': '[#8]=[#6]-[#6]1:[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6]:1-[#6]=[#8]',
        'pms': '[#6]-[#16+](-[#6]-[#6]-[#8]-[#6]=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'esc': '[#8]=[#6]-[#8]-[#6](-[#16](=[#8])(-[#6]-[#6])=[#8])-[#6]',
        'sps': '[#8]=[#6]-[#8]-[#6]-[#6]-[#16](=[#8])(-[#6]1:[#6]:[#6]:[#6](-[#16](=[#8])(-[#8])=[#8]):[#6]:[#6]:1)=[#8]',
    }

    other_alpha_amino_protect_groups_smiles = {
        'benzyloxycarbonyl': 'O=COCC1=CC=CC=C1',
        'allyloxycarbonyl': 'O=COCC=C',
        'o-nitrobenzenesulfonyl': 'O=S(C1=CC=CC=C1[N+]([O-])=O)=O',
        '2,4-dinitrobenzenesulfonyl': 'O=S(C1=CC=C([N+]([O-])=O)C=C1[N+]([O-])=O)=O',
        'benzothiazole-2-sulfonyl': 'O=S(C1=NC2=CC=CC=C2S1)=O',
        '2,2,2-trichloroethyloxycarbonyl': 'O=COCC(Cl)(Cl)Cl',
        'dithiasuccinoyl': 'O=CSS[C]=O',
        'p-nitrobenzyloxycarbonyl': 'O=COCC1=CC=C([N+]([O-])=O)C=C1',
        'alpha-azidoacids': '[N-]=[N+]=NCC(O)=O',
        'proparglyoxycarbonyl': 'C#COC(C)=O',
        'o-nitrobenzylcarbonyl': 'O=CCC1=CC=CC=C1[N+]([O-])=O',
        '4-nitroveratryloxycarbonyl': 'O=COCC1=C([N+]([O-])=O)C=C(OC)C(OC)=C1',
        '2-(2-nitrophenyl)propyloxycarbonyl': 'O=COCC(C1=CC=CC=C1[N+]([O-])=O)C',
        '2-(3,4-methylenedioxy-6-nitrophenyl)propyloxycarbonyl': 'O=COCC(C1=CC(OCO2)=C2C=C1[N+]([O-])=O)C',
        '9-(4-bromophenyl)-9-fluorenyl': 'BrC1=CC=C(C2C3=C(C4=C2C=CC=C4)C=CC=C3)C=C1',
        'azidomethoxycarbonyl': 'O=COCN=[N+]=[N-]',
        'hexafluoroacetone': 'O=C1OC(C(C(F)(F)F)C(F)(F)F)NC1',
    }

    other_alpha_amino_protect_groups_smarts = {
        'benzyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'allyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
        'o-nitrobenzenesulfonyl': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8])=[#8]',
        '2,4-dinitrobenzenesulfonyl': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1-[#7+](-[#8-])=[#8])=[#8]',
        'benzothiazole-2-sulfonyl': '[#8]=[#16](-[#6]1:[#7]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#16]:1)=[#8]',
        '2,2,2-trichloroethyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6](-[#17])(-[#17])-[#17]',
        'dithiasuccinoyl': '[#8]=[#6]-[#16]-[#16]-[#6]=[#8]',
        'p-nitrobenzyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'alpha-azidoacids': '[#7-]=[#7+]=[#7]-[#6]-[#6](-[#8])=[#8]',
        'proparglyoxycarbonyl': '[#6]#[#6]-[#8]-[#6](-[#6])=[#8]',
        'o-nitrobenzylcarbonyl': '[#8]=[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        '4-nitroveratryloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:1',
        '2-(2-nitrophenyl)propyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8])-[#6]',
        '2-(3,4-methylenedioxy-6-nitrophenyl)propyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6](-[#6]1:[#6]:[#6]2-[#8]-[#6]-[#8]-[#6]:2:[#6]:[#6]:1-[#7+](-[#8-])=[#8])-[#6]',
        '9-(4-bromophenyl)-9-fluorenyl': '[#35]-[#6]1:[#6]:[#6]:[#6](-[#6]2-[#6]3:[#6](-[#6]4:[#6]-2:[#6]:[#6]:[#6]:[#6]:4):[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:1',
        'azidomethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#7]=[#7+]=[#7-]',
        'hexafluoroacetone': '[#8]=[#6]1-[#8]-[#6](-[#6](-[#6](-[#9])(-[#9])-[#9])-[#6](-[#9])(-[#9])-[#9])-[#7]-[#6]-1',
    }

    other_alpha_amino_protect_groups_acronym_smiles = {
        'Z': 'O=COCC1=CC=CC=C1',
        'alloc': 'O=COCC=C',
        'o-nbs': 'O=S(C1=CC=CC=C1[N+]([O-])=O)=O',
        'd-nbs': 'O=S(C1=CC=C([N+]([O-])=O)C=C1[N+]([O-])=O)=O',
        'bts': 'O=S(C1=NC2=CC=CC=C2S1)=O',
        'troc': 'O=COCC(Cl)(Cl)Cl',
        'dts': 'O=CSS[C]=O',
        'pnz': 'O=COCC1=CC=C([N+]([O-])=O)C=C1',
        'poc': 'C#COC(C)=O',
        'onz': 'O=CCC1=CC=CC=C1[N+]([O-])=O',
        'nvoc': 'O=COCC1=C([N+]([O-])=O)C=C(OC)C(OC)=C1',
        'nppoc': 'O=COCC(C1=CC=CC=C1[N+]([O-])=O)C',
        'mnppoc': 'O=COCC(C1=CC(OCO2)=C2C=C1[N+]([O-])=O)C',
        'brphf': 'BrC1=CC=C(C2C3=C(C4=C2C=CC=C4)C=CC=C3)C=C1',
        'azoc': 'O=COCN=[N+]=[N-]',
        'hfa': 'O=C1OC(C(C(F)(F)F)C(F)(F)F)NC1',
    }

    other_alpha_amino_protect_groups_acronym_smarts = {
        'Z': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'alloc': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
        'o-nbs': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8])=[#8]',
        'd-nbs': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1-[#7+](-[#8-])=[#8])=[#8]',
        'bts': '[#8]=[#16](-[#6]1:[#7]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#16]:1)=[#8]',
        'troc': '[#8]=[#6]-[#8]-[#6]-[#6](-[#17])(-[#17])-[#17]',
        'dts': '[#8]=[#6]-[#16]-[#16]-[#6]=[#8]',
        'pnz': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'poc': '[#6]#[#6]-[#8]-[#6](-[#6])=[#8]',
        'onz': '[#8]=[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'nvoc': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:1',
        'nppoc': '[#8]=[#6]-[#8]-[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8])-[#6]',
        'mnppoc': '[#8]=[#6]-[#8]-[#6]-[#6](-[#6]1:[#6]:[#6]2-[#8]-[#6]-[#8]-[#6]:2:[#6]:[#6]:1-[#7+](-[#8-])=[#8])-[#6]',
        'brphf': '[#35]-[#6]1:[#6]:[#6]:[#6](-[#6]2-[#6]3:[#6](-[#6]4:[#6]-2:[#6]:[#6]:[#6]:[#6]:4):[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:1',
        'azoc': '[#8]=[#6]-[#8]-[#6]-[#7]=[#7+]=[#7-]',
        'hfa': '[#8]=[#6]1-[#8]-[#6](-[#6](-[#6](-[#9])(-[#9])-[#9])-[#6](-[#9])(-[#9])-[#9])-[#7]-[#6]-1',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_acid_smiles = {
        '2-chlorobenzyloxycarbonyl': 'O=COCC1=CC=CC=C1Cl',
        'tert-butyloxycarbonyl': 'O=COC(C)(C)C',
        '4-methyltrityl': 'CC1=CC=C(C(C2=CC=CC=C2)C3=CC=CC=C3)C=C1',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_acid_smarts = {
        '2-chlorobenzyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#17]',
        'tert-butyloxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        '4-methyltrityl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_acid_acronym_smiles = {
        'cl-z': 'O=COCC1=CC=CC=C1Cl',
        'boc': 'O=COC(C)(C)C',
        'mtt': 'CC1=CC=C(C(C2=CC=CC=C2)C3=CC=CC=C3)C=C1',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_acid_acronym_smarts = {
        'cl-z': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#17]',
        'boc': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        'mtt': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_base_smiles = {
        '9-fluorenylmethoxycarbonyl': 'O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        '1-(4,4-dimethyl-2,6-dioxocylohex-1-ylidene)-3-methylbutyl': r'O=C1/C(C(CC(C)C1)=O)=C\CC(C)C',
        'trifluoroacetyl': 'O=CC(F)(F)F',
        '2-(methylsulfonyl)ethoxycarbonyl': 'O=COCCS(=O)(C)=O',
        'tetrachlorophthaloyl': 'O=CC1=C(Cl)C(Cl)=C(Cl)C(Cl)=C1C=O',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_base_smarts = {
        '9-fluorenylmethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        '1-(4,4-dimethyl-2,6-dioxocylohex-1-ylidene)-3-methylbutyl': '[#8]=[#6]1-[#6](-[#6](-[#6]-[#6](-[#6])-[#6]-1)=[#8])=[#6]-[#6]-[#6](-[#6])-[#6]',
        'trifluoroacetyl': '[#8]=[#6]-[#6](-[#9])(-[#9])-[#9]',
        '2-(methylsulfonyl)ethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]-[#16](=[#8])(-[#6])=[#8]',
        'tetrachlorophthaloyl': '[#8]=[#6]-[#6]1:[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6]:1-[#6]=[#8]',
    }

    lys_orn_dap_dab_protecting_groups_removed_by_base_acronym_smiles = {
        'fmoc': 'O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        'ivdde': r'O=C1/C(C(CC(C)C1)=O)=C\CC(C)C',
        'tfa': 'O=CC(F)(F)F',
        'msc': 'O=COCCS(=O)(C)=O',
        'tcp': 'O=CC1=C(Cl)C(Cl)=C(Cl)C(Cl)=C1C=O'
    }

    lys_orn_dap_dab_protecting_groups_removed_by_base_acronym_smarts = {
        'fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        'ivdde': '[#8]=[#6]1-[#6](-[#6](-[#6]-[#6](-[#6])-[#6]-1)=[#8])=[#6]-[#6]-[#6](-[#6])-[#6]',
        'tfa': '[#8]=[#6]-[#6](-[#9])(-[#9])-[#9]',
        'msc': '[#8]=[#6]-[#8]-[#6]-[#6]-[#16](=[#8])(-[#6])=[#8]',
        'tcp': '[#8]=[#6]-[#6]1:[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6](-[#17]):[#6]:1-[#6]=[#8]',
    }

    other_lys_orn_dap_dab_protecting_groups_smiles = {
        'allyloxycarbonyl': 'O=COCC=C',
        'p-nitrobenzyloxycarbonyl': 'O=COCC1=CC=C([N+]([O-])=O)C=C1',
        'phenyldisulphanylethyloxycarbonyl': 'O=COC(SSC1=CC=CC=C1)C',
        '2-pyridyldisulphanylethyloxycarbonyl': 'O=COC(SSC1=NC=CC=C1)C',
        'o-nitrobenzenesulfonyl': 'O=S(C1=CC=CC=C1[N+]([O-])=O)=O',
    }

    other_lys_orn_dap_dab_protecting_groups_smarts = {
        'allyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
        'p-nitrobenzyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'phenyldisulphanylethyloxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#16]-[#16]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]',
        '2-pyridyldisulphanylethyloxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#16]-[#16]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1)-[#6]',
        'o-nitrobenzenesulfonyl': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8])=[#8]',
    }

    other_lys_orn_dap_dab_protecting_groups_acronym_smiles = {
        'alloc': 'O=COCC=C',
        'pnz': 'O=COCC1=CC=C([N+]([O-])=O)C=C1',
        'phdec': 'O=COC(SSC1=CC=CC=C1)C',
        'pydec': 'O=COC(SSC1=NC=CC=C1)C',
        'o-nbs': 'O=S(C1=CC=CC=C1[N+]([O-])=O)=O'
    }

    other_lys_orn_dap_dab_protecting_groups_acronym_smarts = {
        'alloc': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
        'pnz': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'phdec': '[#8]=[#6]-[#8]-[#6](-[#16]-[#16]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]',
        'pydec': '[#8]=[#6]-[#8]-[#6](-[#16]-[#16]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1)-[#6]',
        'o-nbs': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8])=[#8]',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_acid_smiles = {
        'tert-butyl': 'CC(C)C',
        '2-chlorotrityl': 'ClC1=CC=CC=C1C(C2=CC=CC=C2)C3=CC=CC=C3',
        '2-4-dimethyoxybenzyl': 'CC1=CC(OC)=CC(OC)=C1',
        '2-phenylisopropyl': 'CC(C)C1=CC=CC=C1',
        '5-phenyl-3,4-ethylenedioxythenyl': 'CC1=C(OCCO2)C2=C(S1)C3=CC=CC=C3',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_acid_smarts = {
        'tert-butyl': '[#6]-[#6](-[#6])-[#6]',
        '2-chlorotrityl': '[#17]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '2-4-dimethyoxybenzyl': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:1',
        '2-phenylisopropyl': '[#6]-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '5-phenyl-3,4-ethylenedioxythenyl': '[#6]-[#6]1:[#6]2-[#8]-[#6]-[#6]-[#8]-[#6]:2:[#6](:[#16]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_acid_acronym_smiles = {
        'bu': 'CC(C)C',
        '2-cl-trt': 'ClC1=CC=CC=C1C(C2=CC=CC=C2)C3=CC=CC=C3',
        'dmb': 'CC1=CC(OC)=CC(OC)=C1',
        '2-ph-pr': 'CC(C)C1=CC=CC=C1',
        'phenyl-edotn': 'CC1=C(OCCO2)C2=C(S1)C3=CC=CC=C3',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_acid_acronym_smarts = {
        'bu': '[#6]-[#6](-[#6])-[#6]',
        '2-cl-trt': '[#17]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'dmb': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:1',
        '2-ph-pr': '[#6]-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'phenyl-edotn': '[#6]-[#6]1:[#6]2-[#8]-[#6]-[#6]-[#8]-[#6]:2:[#6](:[#16]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_base_smiles = {
        '9-fluorenylmethyl': 'CC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        '4-(N-[1-(4,4-dimethyl-2,6-dioxocylocheylidene)-3-methylbutyl]-amino)benzyl': r'CC(CC(/C1=C(NC2=CC=C(C)C=C2)/CC(C)C)=O)(C)CC1=O',
        'methyl': 'C',
        'ethyl': 'CC',
        'carbamoylmethyl': 'CC(N)=O',

    }

    alpha_carboxylic_acid_protecting_groups_removed_by_base_smarts = {
        '9-fluorenylmethyl': '[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        '4-(N-[1-(4,4-dimethyl-2,6-dioxocylocheylidene)-3-methylbutyl]-amino)benzyl': '[#6]-[#6]1(-[#6]-[#6](-[#6](=[#6](-[#7]-[#6]2:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:2)-[#6]-[#6](-[#6])-[#6])-[#6](-[#6]-1)=[#8])=[#8])-[#6]',
        'methyl': '[#6]',
        'ethyl': '[#6]-[#6]',
        'carbamoylmethyl': '[#6]-[#6](-[#7])=[#8]',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_base_acronym_smiles = {
        'fm': 'CC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        'dmab': r'CC(CC(/C1=C(NC2=CC=C(C)C=C2)/CC(C)C)=O)(C)CC1=O',
        'me': 'C',
        'et': 'CC',
        'cam': 'CC(N)=O',
    }

    alpha_carboxylic_acid_protecting_groups_removed_by_base_acronym_smarts = {
        'fm': '[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        'dmab': '[#6]-[#6]1(-[#6]-[#6](-[#6](=[#6](-[#7]-[#6]2:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:2)-[#6]-[#6](-[#6])-[#6])-[#6](-[#6]-1)=[#8])=[#8])-[#6]',
        'me': '[#6]',
        'et': '[#6]-[#6]',
        'cam': '[#6]-[#6](-[#7])=[#8]',
    }

    other_alpha_carboxylic_acid_protecting_group_smiles = {
        'allyl': 'CC=C',
        'benzyl': 'CC1=CC=CC=C1',
        'phenacyl': 'CC(C1=CC=CC=C1)=O',
        'p-nitrobenzyl': 'CC1=CC=C([N+]([O-])=O)C=C1',
        '2-trimethylsilyethyl': 'CCC[Si](C)(C)C',
        '(2-phenyl-2-trimethylsilyl)ethyl': 'CC(C1=CC=CC=C1)[Si](C)(C)C',
        '2-(trimethylsilyl)isopropyl': 'CC(C)([Si](C)(C)C)C',
        '2,2,2-trichloroethyl': 'CC(Cl)(Cl)Cl',
        'p-hydroxyphenacyl': 'CC(C1=CC=C(O)C=C1)=O',
        '4,5-dimethyoxy-2-nitrobenzyl': 'CC1=CC(OC)=C(OC)C=C1[N+]([O-])=O',
        '1,1-dimethylallyl': 'C=CC(C)C',
        'pentaaminecobalt_III': 'N[Co](N)(N)(N)(Cl)(Cl)N',
    }

    other_alpha_carboxylic_acid_protecting_group_smarts = {
        'allyl': '[#6]-[#6]=[#6]',
        'benzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'phenacyl': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8]',
        'p-nitrobenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        '2-trimethylsilyethyl': '[#6]-[#6]-[#6]-[Si](-[#6])(-[#6])-[#6]',
        '(2-phenyl-2-trimethylsilyl)ethyl': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[Si](-[#6])(-[#6])-[#6]',
        '2-(trimethylsilyl)isopropyl': '[#6]-[#6](-[#6])(-[Si](-[#6])(-[#6])-[#6])-[#6]',
        '2,2,2-trichloroethyl': '[#6]-[#6](-[#17])(-[#17])-[#17]',
        'p-hydroxyphenacyl': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6](-[#8]):[#6]:[#6]:1)=[#8]',
        '4,5-dimethyoxy-2-nitrobenzyl': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        '1,1-dimethylallyl': '[#6]=[#6]-[#6](-[#6])-[#6]',
        'pentaaminecobalt_III': '[#7]-[Co](-[#7])(-[#7])(-[#7])(-[#17])(-[#17])-[#7]',
    }

    other_alpha_carboxylic_acid_protecting_group_acronym_smiles = {
        'al': 'CC=C',
        'bn': 'CC1=CC=CC=C1',
        'pac': 'CC(C1=CC=CC=C1)=O',
        'pnb': 'CC1=CC=C([N+]([O-])=O)C=C1',
        'tmse': 'CCC[Si](C)(C)C',
        'ptmse': 'CC(C1=CC=CC=C1)[Si](C)(C)C',
        'tmsi': 'CC(C)([Si](C)(C)C)C',
        'tce': 'CC(Cl)(Cl)Cl',
        'php': 'CC(C1=CC=C(O)C=C1)=O',
        'dmnb': 'CC1=CC(OC)=C(OC)C=C1[N+]([O-])=O',
        'dma': 'C=CC(C)C',
    }

    other_alpha_carboxylic_acid_protecting_group_acronym_smarts = {
        'al': '[#6]-[#6]=[#6]',
        'bn': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'pac': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8]',
        'pnb': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'tmse': '[#6]-[#6]-[#6]-[Si](-[#6])(-[#6])-[#6]',
        'ptmse': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[Si](-[#6])(-[#6])-[#6]',
        'tmsi': '[#6]-[#6](-[#6])(-[Si](-[#6])(-[#6])-[#6])-[#6]',
        'tce': '[#6]-[#6](-[#17])(-[#17])-[#17]',
        'php': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6](-[#8]):[#6]:[#6]:1)=[#8]',
        'dmnb': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'dma': '[#6]=[#6]-[#6](-[#6])-[#6]',
    }

    asp_glu_protecting_groups_removed_by_acid_smiles = {
        'benzyl': 'CC1=CC=CC=C1',
        'cyclohexyl': 'C1CCCCC1',
        'tert-butyl': 'CC(C)C',
        'b-menthyl': 'C[C@H]1[C@H](C(C)C)CC[C@@H](C)C1',
        'b-3-methylpent-3-yl': 'CCC(C)CC',
        '2-phenylisopropyl': 'CC(C)C1=CC=CC=C1',
        '4-(3,6,9-trioxadecyl)oxybenzyl': 'CC1=CC=C(OCCOCCOCCOC)C=C1',
    }

    asp_glu_protecting_groups_removed_by_acid_smarts = {
        'benzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'cyclohexyl': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        'tert-butyl': '[#6]-[#6](-[#6])-[#6]',
        'b-menthyl': '[#6]-[#6@H]1-[#6@H](-[#6](-[#6])-[#6])-[#6]-[#6]-[#6@@H](-[#6])-[#6]-1',
        'b-3-methylpent-3-yl': '[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
        '2-phenylisopropyl': '[#6]-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '4-(3,6,9-trioxadecyl)oxybenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]):[#6]:[#6]:1',
    }

    asp_glu_protecting_groups_removed_by_acid_acronym_smiles = {
        'bn': 'CC1=CC=CC=C1',
        'chx': 'C1CCCCC1',
        'bu': 'CC(C)C',
        'men': 'C[CH]1[CH](C(C)C)CC[CH](C)C1',
        'mpe': 'CCC(C)CC',
        '2-ph-pr': 'CC(C)C1=CC=CC=C1',
        'tegbz': 'CC1=CC=C(OCCOCCOCCOC)C=C1',
    }

    asp_glu_protecting_groups_removed_by_acid_acronym_smarts = {
        'bn': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'chx': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        'bu': '[#6]-[#6](-[#6])-[#6]',
        'men': '[#6]-[#6H]1-[#6H](-[#6](-[#6])-[#6])-[#6]-[#6]-[#6H](-[#6])-[#6]-1',
        'mpe': '[#6]-[#6]-[#6](-[#6])-[#6]-[#6]',
        '2-ph-pr': '[#6]-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tegbz': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]):[#6]:[#6]:1',
    }

    asp_glu_protecting_groups_removed_by_base_smiles = {
        '9-fluoroenylmethyl': 'CC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        '4-(N-[1-(4,4-dimethyl-2,6-dioxocyclohexylidene)-3-methyl-butyl]-amino)benzyl': r'CC(CC(/C1=C(NC2=CC=C(C)C=C2)/CC(C)C)=O)(C)CC1=O',
    }

    asp_glu_protecting_groups_removed_by_base_smarts = {
        '9-fluoroenylmethyl': '[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        '4-(N-[1-(4,4-dimethyl-2,6-dioxocyclohexylidene)-3-methyl-butyl]-amino)benzyl': '[#6]-[#6]1(-[#6]-[#6](-[#6](=[#6](-[#7]-[#6]2:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:2)-[#6]-[#6](-[#6])-[#6])-[#6](-[#6]-1)=[#8])=[#8])-[#6]',
    }

    asp_glu_protecting_groups_removed_by_base_acronym_smiles = {
        'fm': 'CC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        'dmab': r'CC(CC(/C1=C(NC2=CC=C(C)C=C2)/CC(C)C)=O)(C)CC1=O'
    }

    asp_glu_protecting_groups_removed_by_base_acronym_smarts = {
        'fm': '[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        'dmab': '[#6]-[#6]1(-[#6]-[#6](-[#6](=[#6](-[#7]-[#6]2:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:2)-[#6]-[#6](-[#6])-[#6])-[#6](-[#6]-1)=[#8])=[#8])-[#6]',
    }

    other_asp_glu_protecting_groups_smiles = {
        'allyl': 'CC=C',
        'p-nitrobenzyl': 'CC1=CC=C([N+]([O-])=O)C=C1',
        'trimethylsilylethyl': 'CCC[Si](C)(C)C',
        '(2-phenyl-2-trimethylsilyl)ethyl': 'CC(C1=CC=CC=C1)[Si](C)(C)C',
        '4,5-dimethoxy-2-nitrobenzyloxycarbonyl': 'CC1=CC(OC)=C(OC)C=C1[N+]([O-])=O',
    }

    other_asp_glu_protecting_groups_smarts = {
        'allyl': '[#6]-[#6]=[#6]',
        'p-nitrobenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'trimethylsilylethyl': '[#6]-[#6]-[#6]-[Si](-[#6])(-[#6])-[#6]',
        '(2-phenyl-2-trimethylsilyl)ethyl': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[Si](-[#6])(-[#6])-[#6]',
        '4,5-dimethoxy-2-nitrobenzyloxycarbonyl': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
    }

    other_asp_glu_protecting_groups_acronym_smiles = {
        'al': 'CC=C',
        'pnb': 'CC1=CC=C([N+]([O-])=O)C=C1',
        'tmse': 'CCC[Si](C)(C)C',
        'ptmse': 'CC(C1=CC=CC=C1)[Si](C)(C)C',
        'dmnb': 'CC1=CC(OC)=C(OC)C=C1[N+]([O-])=O',
    }

    other_asp_glu_protecting_groups_acronym_smarts = {
        'al': '[#6]-[#6]=[#6]',
        'pnb': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1',
        'tmse': '[#6]-[#6]-[#6]-[Si](-[#6])(-[#6])-[#6]',
        'ptmse': '[#6]-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[Si](-[#6])(-[#6])-[#6]',
        'dmnb': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
    }

    amide_backbone_protecting_group_removed_by_acid_smiles = {
        'pseudoprolines': 'CC1(C)NC(C(O)=O)CO1',
        '2-hydroxy-4-methoxybenzyl': 'CC1=CC=C(OC)C=C1O',
        '2,4-dimethoxybenzyl': 'CC1=CC=C(OC)C=C1OC',
        '2,4,6-trimethoxybenzyl': 'CC1=C(OC)C=C(OC)C=C1OC',
        '1-methyl-3-indolylmethyl': 'CCC1=CNC2=C1C=CC=C2',
        '3,4-ethylene-dioxy-2-thenyl': 'CC1=C(OCCO2)C2=CS1',
    }

    amide_backbone_protecting_group_removed_by_acid_smarts = {
        'pseudoprolines': '[#6]-[#6]1(-[#6])-[#7]-[#6](-[#6](-[#8])=[#8])-[#6]-[#8]-1',
        '2-hydroxy-4-methoxybenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]',
        '2,4-dimethoxybenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
        '2,4,6-trimethoxybenzyl': '[#6]-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
        '1-methyl-3-indolylmethyl': '[#6]-[#6]-[#6]1:[#6]:[#7H]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
        '3,4-ethylene-dioxy-2-thenyl': '[#6]-[#6]1:[#6]2-[#8]-[#6]-[#6]-[#8]-[#6]:2:[#6]:[#16]:1',
    }

    amide_backbone_protecting_group_removed_by_acid_acronym_smiles = {
        'hmb': 'CC1=CC=C(OC)C=C1O',
        'dmb': 'CC1=CC=C(OC)C=C1OC',
        'tmob': 'CC1=C(OC)C=C(OC)C=C1OC',
        'mim': 'CCC1=CNC2=C1C=CC=C2',
        'edot': 'CC1=C(OCCO2)C2=CS1',
    }

    amide_backbone_protecting_group_removed_by_acid_acronym_smarts  = {
        'hmb': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]',
        'dmb': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
        'tmob': '[#6]-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
        'mim': '[#6]-[#6]-[#6]1:[#6]:[#7H]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
        'edot': '[#6]-[#6]1:[#6]2-[#8]-[#6]-[#6]-[#8]-[#6]:2:[#6]:[#16]:1',
    }

    other_amide_backbone_protecting_group_smiles = {
        '4-methoxy-2-nitro-benzyl': 'CC1=CC=C(OC)C=C1[N+]([O-])=O',
        '(6-hydroxy-3-oxido-1,3-benz[d]oxathiol-5-yl)methyl': 'O=S1COC2=C1C=C(C)C(O)=C2',
        '2-hydroxy-4-methoxy-5-(methylsulfinyl)benzyl': 'CC1=CC(S(C)=O)=C(OC)C=C1O',
        'n-boc-n-methyl[2-(methylamino)ethyl]carbamoyl-hmb': 'CC(C)(OC(N(CCN(C(OC1=CC(OC)=CC=C1)=O)C)C)=O)C',
    }

    other_amide_backbone_protecting_group_smarts = {
        '4-methoxy-2-nitro-benzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        '(6-hydroxy-3-oxido-1,3-benz[d]oxathiol-5-yl)methyl': '[#8]=[#16]1-[#6]-[#8]-[#6]2:[#6]-1:[#6]:[#6](-[#6]):[#6](-[#8]):[#6]:2',
        '2-hydroxy-4-methoxy-5-(methylsulfinyl)benzyl': '[#6]-[#6]1:[#6]:[#6](-[#16](-[#6])=[#8]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]',
        'n-boc-n-methyl[2-(methylamino)ethyl]carbamoyl-hmb': '[#6]-[#6](-[#6])(-[#8]-[#6](-[#7](-[#6]-[#6]-[#7](-[#6](-[#8]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:[#6]:1)=[#8])-[#6])-[#6])=[#8])-[#6]',
    }

    asn_gln_protecting_groups_removed_by_acid_smiles = {
        '9-xanthenyl': 'C1(CC2=C(C=CC=C2)O3)=C3C=CC=C1',
        'trityl': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        '4-methyltrityl': 'CC1=CC=C(C(C2=CC=CC=C2)C3=CC=CC=C3)C=C1',
        'cyclopropyldimethylcarbinyl': 'CC(C1CC1)C',
        '4,4-dimethoxybenzhydryl': 'COC1=CC=C(C=C1)CC2=CC=C(OC)C=C2',
        '2,4,6-trimethoxybenzyl': 'CC1=C(OC)C=C(OC)C=C1OC',
    }

    asn_gln_protecting_groups_removed_by_acid_smarts = {
        '9-xanthenyl': '[#6]12-[#6]-[#6]3:[#6](:[#6]:[#6]:[#6]:[#6]:3)-[#8]-[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
        'trityl': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '4-methyltrityl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
        'cyclopropyldimethylcarbinyl': '[#6]-[#6](-[#6]1-[#6]-[#6]-1)-[#6]',
        '4,4-dimethoxybenzhydryl': '[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1',
        '2,4,6-trimethoxybenzyl': '[#6]-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
    }

    asn_gln_protecting_groups_removed_by_acid_acronym_smiles = {
        'xan': 'C1(CC2=C(C=CC=C2)O3)=C3C=CC=C1',
        'trt': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        'mtt': 'CC1=CC=C(C(C2=CC=CC=C2)C3=CC=CC=C3)C=C1',
        'cpd': 'CC(C1CC1)C',
        'mbh': 'COC1=CC=C(C=C1)CC2=CC=C(OC)C=C2',
        'tmob': 'CC1=C(OC)C=C(OC)C=C1OC',
    }

    asn_gln_protecting_groups_removed_by_acid_acronym_smarts = {
        'xan': '[#6]12-[#6]-[#6]3:[#6](:[#6]:[#6]:[#6]:[#6]:3)-[#8]-[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
        'trt': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'mtt': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
        'cpd': '[#6]-[#6](-[#6]1-[#6]-[#6]-1)-[#6]',
        'mbh': '[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1',
        'tmob': '[#6]-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
    }

    arg_protecting_groups_removed_by_acid_smiles = {
        'p-toluenesulfonyl': 'O=S(C1=CC=C(C)C=C1)=O',
        '2,2,5,7,8-pentamethylchroman-6-sulfonyl': 'O=S(C1=C(C)C(CCC(C)(C)O2)=C2C(C)=C1C)=O',
        '2,2,4,6,7-pentamethyl-2,3-dihydrobenzofuran-5-sulfonyl': 'O=S(C1=C(C)C(C)=C(OC(C)(C)C2)C2=C1C)=O',
        'mesityl-2-sulfonyl': 'CC1=C(S(=O)(NC(N)=N)=O)C(C)=CC(C)=C1',
        '4-methoxy-2,3,6-trimethylphenylsulfonyl': 'O=S(C1=C(C)C=C(OC)C(C)=C1C)=O',
        '1,2-dimethylindole-3-sulfonyl': 'O=S(C1=C(C)N(C)C2=C1C=CC=C2)=O',
        'w,w-bis-tert-butyloxycarbonyl': r'CC(C)(OC(/N=C(NC(OC(C)(C)C)=O)\N)=O)C',
        '5-dibenzosuberenyl': 'C12=CC=CC=C1CCC3C=CC=CC3=C2',
        '5-dibenzosuberyl': 'C12=CC=CC=C1CCC3=C(C=CC=C3)C2',
        '2-methoxy-5-dibenzosuberyl': 'COC(C=C1CC2)=CC=C1CC3=C2C=CC=C3',
        'nitro': 'O=[N+][O-]'
    }

    arg_protecting_groups_removed_by_acid_smarts = {
        'p-toluenesulfonyl': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1)=[#8]',
        '2,2,5,7,8-pentamethylchroman-6-sulfonyl': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#6]2-[#6]-[#6]-[#6](-[#6])(-[#6])-[#8]-[#6]:2:[#6](-[#6]):[#6]:1-[#6])=[#8]',
        '2,2,4,6,7-pentamethyl-2,3-dihydrobenzofuran-5-sulfonyl': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#6](-[#6]):[#6]2-[#8]-[#6](-[#6])(-[#6])-[#6]-[#6]:2:[#6]:1-[#6])=[#8]',
        'mesityl-2-sulfonyl': '[#6]-[#6]1:[#6](-[#16](=[#8])(-[#7]-[#6](-[#7])=[#7])=[#8]):[#6](-[#6]):[#6]:[#6](-[#6]):[#6]:1',
        '4-methoxy-2,3,6-trimethylphenylsulfonyl': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#6]:[#6](-[#8]-[#6]):[#6](-[#6]):[#6]:1-[#6])=[#8]',
        '1,2-dimethylindole-3-sulfonyl': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#7](-[#6]):[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2)=[#8]',
        'w,w-bis-tert-butyloxycarbonyl': '[#6]-[#6](-[#6])(-[#8]-[#6](/[#7]=[#6](/[#7]-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])=[#8])-[#7])=[#8])-[#6]',
        '5-dibenzosuberenyl': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]1-[#6]=[#6]-[#6]=[#6]-[#6]-1=[#6]-2',
        '5-dibenzosuberyl': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]1:[#6](:[#6]:[#6]:[#6]:[#6]:1)-[#6]-2',
        '2-methoxy-5-dibenzosuberyl': '[#6]-[#8]-[#6]1:[#6]:[#6]2-[#6]-[#6]-[#6]3:[#6](-[#6]-[#6]:2:[#6]:[#6]:1):[#6]:[#6]:[#6]:[#6]:3',
        'nitro': '[#8]=[#7+]-[#8-]',
    }

    arg_protecting_groups_removed_by_acid_acronym_smiles = {
        'tos': 'O=S(C1=CC=C(C)C=C1)=O',
        'pmc': 'O=S(C1=C(C)C(CCC(C)(C)O2)=C2C(C)=C1C)=O',
        'pbf': 'O=S(C1=C(C)C(C)=C(OC(C)(C)C2)C2=C1C)=O',
        'mts': 'CC1=C(S(=O)(NC(N)=N)=O)C(C)=CC(C)=C1',
        'mtr': 'O=S(C1=C(C)C=C(OC)C(C)=C1C)=O',
        'mis': 'O=S(C1=C(C)N(C)C2=C1C=CC=C2)=O',
        'bis-boc': r'CC(C)(OC(/N=C(NC(OC(C)(C)C)=O)\N)=O)C',
        'suben': 'C12=CC=CC=C1CCC3C=CC=CC3=C2',
        'sub': 'C12=CC=CC=C1CCC3=C(C=CC=C3)C2',
        'mesub': 'COC(C=C1CC2)=CC=C1CC3=C2C=CC=C3',
        'no2': 'O=[N+][O-]'
    }

    arg_protecting_groups_removed_by_acid_acronym_smarts = {
        'tos': '[#8]=[#16](-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1)=[#8]',
        'pmc': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#6]2-[#6]-[#6]-[#6](-[#6])(-[#6])-[#8]-[#6]:2:[#6](-[#6]):[#6]:1-[#6])=[#8]',
        'pbf': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#6](-[#6]):[#6]2-[#8]-[#6](-[#6])(-[#6])-[#6]-[#6]:2:[#6]:1-[#6])=[#8]',
        'mts': '[#6]-[#6]1:[#6](-[#16](=[#8])(-[#7]-[#6](-[#7])=[#7])=[#8]):[#6](-[#6]):[#6]:[#6](-[#6]):[#6]:1',
        'mtr': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#6]:[#6](-[#8]-[#6]):[#6](-[#6]):[#6]:1-[#6])=[#8]',
        'mis': '[#8]=[#16](-[#6]1:[#6](-[#6]):[#7](-[#6]):[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2)=[#8]',
        'bis-boc': '[#6]-[#6](-[#6])(-[#8]-[#6](/[#7]=[#6](/[#7]-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])=[#8])-[#7])=[#8])-[#6]',
        'suben': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]1-[#6]=[#6]-[#6]=[#6]-[#6]-1=[#6]-2',
        'sub': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]1:[#6](:[#6]:[#6]:[#6]:[#6]:1)-[#6]-2',
        'mesub': '[#6]-[#8]-[#6]1:[#6]:[#6]2-[#6]-[#6]-[#6]3:[#6](-[#6]-[#6]:2:[#6]:[#6]:1):[#6]:[#6]:[#6]:[#6]:3',
        'no2': '[#8]=[#7+]-[#8-]',
    }

    arg_protecting_groups_removed_by_base_smiles = {
        'trifluoroacetyl': 'O=CC(F)(F)F'
    }

    arg_protecting_groups_removed_by_base_smarts = {
        'trifluoroacetyl': '[#8]=[#6]-[#6](-[#9])(-[#9])-[#9]',
    }

    arg_protecting_groups_removed_by_base_acronym_smiles = {
        'tfa': 'O=CC(F)(F)F'
    }

    arg_protecting_groups_removed_by_base_acronym_smarts = {
        'tfa': '[#8]=[#6]-[#6](-[#9])(-[#9])-[#9]',
    }

    other_arg_protecting_groups_smiles = {
        'w,w-bis-benzyloxycarbonyl': r'O=C(/N=C(NC(OCC1=CC=CC=C1)=O)\N)OCC2=CC=CC=C2',
        'w,w-bis-allyloxycarbonyl': r'O=C(/N=C(NC(OCC=C)=O)\N)OCC=C'
    }

    other_arg_protecting_groups_smarts = {
        'w,w-bis-benzyloxycarbonyl': '[#8]=[#6](/[#7]=[#6](/[#7]-[#6](-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8])-[#7])-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'w,w-bis-allyloxycarbonyl': '[#8]=[#6](/[#7]=[#6](/[#7]-[#6](-[#8]-[#6]-[#6]=[#6])=[#8])-[#7])-[#8]-[#6]-[#6]=[#6]',
    }

    other_arg_protecting_groups_acronym_smiles = {
        'z': r'O=C(/N=C(NC(OCC1=CC=CC=C1)=O)\N)OCC2=CC=CC=C2',
        'alloc': r'O=C(/N=C(NC(OCC=C)=O)\N)OCC=C',
    }

    other_arg_protecting_groups_acronym_smarts = {
        'z': '[#8]=[#6](/[#7]=[#6](/[#7]-[#6](-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8])-[#7])-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'alloc': '[#8]=[#6](/[#7]=[#6](/[#7]-[#6](-[#8]-[#6]-[#6]=[#6])=[#8])-[#7])-[#8]-[#6]-[#6]=[#6]',
    }

    cys_protecting_groups_removed_by_acid_smiles = {
        'p-methylbenzyl': 'CC1=CC=C(C)C=C1',
        'p-methoxybenzyl': 'CC1=CC=C(OC)C=C1',
        'trityl': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        'monomethoxytrityl': 'COC1=CC=C(C(C2=CC=CC=C2)C3=CC=CC=C3)C=C1',
        'trimethoxybenzyl': 'CC1=C(OC)C=C(OC)C=C1OC',
        '9-xanthenyl': 'C1(CC2=C(C=CC=C2)O3)=C3C=CC=C1',
        '2,2,4,6,7-pentamethyl-5-dihydrobenzofuranylmethyl': 'CC1=C(C)C(C)=C(OC(C)(C)C2)C2=C1C',
        'benzyl': 'CC1=CC=CC=C1',
        'tert-butyl': 'CC(C)C',
        '1-adamantyl': 'C12CC3CC(C2)CC(C3)C1',
    }

    cys_protecting_groups_removed_by_acid_smarts = {
        'p-methylbenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1',
        'p-methoxybenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1',
        'trityl': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'monomethoxytrityl': '[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
        'trimethoxybenzyl': '[#6]-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
        '9-xanthenyl': '[#6]12-[#6]-[#6]3:[#6](:[#6]:[#6]:[#6]:[#6]:3)-[#8]-[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
        '2,2,4,6,7-pentamethyl-5-dihydrobenzofuranylmethyl': '[#6]-[#6]1:[#6](-[#6]):[#6](-[#6]):[#6]2-[#8]-[#6](-[#6])(-[#6])-[#6]-[#6]:2:[#6]:1-[#6]',
        'benzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tert-butyl': '[#6]-[#6](-[#6])-[#6]',
        '1-adamantyl': '[#6]12-[#6]-[#6]3-[#6]-[#6](-[#6]-1)-[#6]-[#6](-[#6]-3)-[#6]-2',
    }

    cys_protecting_groups_removed_by_acid_acronym_smiles = {
        'meb': 'CC1=CC=C(C)C=C1',
        'mob': 'CC1=CC=C(OC)C=C1',
        'trt': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        'mmt': 'COC1=CC=C(C(C2=CC=CC=C2)C3=CC=CC=C3)C=C1',
        'tmob': 'CC1=C(OC)C=C(OC)C=C1OC',
        'xan': 'C1(CC2=C(C=CC=C2)O3)=C3C=CC=C1',
        'pmbf': 'CC1=C(C)C(C)=C(OC(C)(C)C2)C2=C1C',
        'bn': 'CC1=CC=CC=C1',
        'bu': 'CC(C)C',
        '1-ada': 'C12CC3CC(C2)CC(C3)C1',
    }

    cys_protecting_groups_removed_by_acid_acronym_smarts = {
        'meb': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1',
        'mob': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1',
        'trt': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'mmt': '[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:1',
        'tmob': '[#6]-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6](-[#8]-[#6]):[#6]:[#6]:1-[#8]-[#6]',
        'xan': '[#6]12-[#6]-[#6]3:[#6](:[#6]:[#6]:[#6]:[#6]:3)-[#8]-[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
        'pmbf': '[#6]-[#6]1:[#6](-[#6]):[#6](-[#6]):[#6]2-[#8]-[#6](-[#6])(-[#6])-[#6]-[#6]:2:[#6]:1-[#6]',
        'bn': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'bu': '[#6]-[#6](-[#6])-[#6]',
        '1-ada': '[#6]12-[#6]-[#6]3-[#6]-[#6](-[#6]-1)-[#6]-[#6](-[#6]-3)-[#6]-2',
    }

    cys_protecting_groups_removed_by_base_smiles = {
        '9-fluorenylmethyl': 'CC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        '2-(2,4-dinitrophenyl)ethyl': 'CCC1=CC=C([N+]([O-])=O)C=C1[N+]([O-])=O',
        '9-fluororenylmethoxycarbonyl': 'O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2',
    }

    cys_protecting_groups_removed_by_base_smarts = {
        '9-fluorenylmethyl': '[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        '2-(2,4-dinitrophenyl)ethyl': '[#6]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        '9-fluororenylmethoxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
    }

    cys_protecting_groups_removed_by_base_acronym_smiles = {
        'fm': 'CC1C2=C(C3=C1C=CC=C3)C=CC=C2',
        'dnpe': 'CCC1=CC=C([N+]([O-])=O)C=C1[N+]([O-])=O',
        'fmoc': 'O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2',
    }

    cys_protecting_groups_removed_by_base_acronym_smarts = {
        'fm': '[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
        'dnpe': '[#6]-[#6]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'fmoc': '[#8]=[#6]-[#8]-[#6]-[#6]1-[#6]2:[#6](-[#6]3:[#6]-1:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2',
    }

    other_cys_protecting_groups_smiles = {
        'acetamidomethyl': 'CNC(C)=O',
        'phenylacetamidomethyl': 'CNC(CC1=CC=CC=C1)=O',
        '5-tert-butylmercapto': 'CC(C)(S)C',
        '3-nitro-2-pyridinesulfenyl': 'SC1=NC=CC=C1[N+]([O-])=O',
        '2-pyridinesulfenyl': 'SC1=NC=CC=C1',
        'allyloxycarbonyl': 'O=COCC=C',
        'N-allyloxycarbonyl-N-[2,3,5,6-tetrafluoro-4-(phenylthio)phenyl]]aminomethyl': 'FC1=C(F)C(SC2=CC=CC=C2)=C(F)C(F)=C1N(C(OCC=C)=O)C',
        'o-nitrobenzyl': 'CC1=CC=CC=C1[N+]([O-])=O',
        '4-picolyl': 'CC1=CC=NC=C1',
        'ninhydrin': 'O=C1C2(SCC(C(O)=O)N2)C(C3=C1C=CC=C3)=O',
    }

    other_cys_protecting_groups_smarts = {
        'acetamidomethyl': '[#6]-[#7]-[#6](-[#6])=[#8]',
        'phenylacetamidomethyl': '[#6]-[#7]-[#6](-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8]',
        '5-tert-butylmercapto': '[#6]-[#6](-[#6])(-[#16])-[#6]',
        '3-nitro-2-pyridinesulfenyl': '[#16]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        '2-pyridinesulfenyl': '[#16]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1',
        'allyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
        'N-allyloxycarbonyl-N-[2,3,5,6-tetrafluoro-4-(phenylthio)phenyl]]aminomethyl': '[#9]-[#6]1:[#6](-[#9]):[#6](-[#16]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6](-[#9]):[#6](-[#9]):[#6]:1-[#7](-[#6](-[#8]-[#6]-[#6]=[#6])=[#8])-[#6]',
        'o-nitrobenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        '4-picolyl': '[#6]-[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1',
        'ninhydrin': '[#8]=[#6]1-[#6]2(-[#16]-[#6]-[#6](-[#6](-[#8])=[#8])-[#7]-2)-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#8]',
    }

    other_cys_protecting_groups_acronym_smiles = {
        'acm': 'CNC(C)=O',
        'phacm': 'CNC(CC1=CC=CC=C1)=O',
        'sbu': 'CC(C)(S)C',
        'npys': 'SC1=NC=CC=C1[N+]([O-])=O',
        's-pyr': 'SC1=NC=CC=C1',
        'alloc': 'O=COCC=C',
        'fsam': 'FC1=C(F)C(SC2=CC=CC=C2)=C(F)C(F)=C1N(C(OCC=C)=O)C',
        'onb': 'CC1=CC=CC=C1[N+]([O-])=O',
        'nin': 'O=C1C2(SCC(C(O)=O)N2)C(C3=C1C=CC=C3)=O'
    }

    other_cys_protecting_groups_acronym_smarts = {
        'acm': '[#6]-[#7]-[#6](-[#6])=[#8]',
        'phacm': '[#6]-[#7]-[#6](-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)=[#8]',
        'sbu': '[#6]-[#6](-[#6])(-[#16])-[#6]',
        'npys': '[#16]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        's-pyr': '[#16]-[#6]1:[#7]:[#6]:[#6]:[#6]:[#6]:1',
        'alloc': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
        'fsam': '[#9]-[#6]1:[#6](-[#9]):[#6](-[#16]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6](-[#9]):[#6](-[#9]):[#6]:1-[#7](-[#6](-[#8]-[#6]-[#6]=[#6])=[#8])-[#6]',
        'onb': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'nin': '[#8]=[#6]1-[#6]2(-[#16]-[#6]-[#6](-[#6](-[#8])=[#8])-[#7]-2)-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#8]',
    }

    his_protecting_groups_removed_by_acid_smiles = {
        'n-tosyl': 'O=S(N1C=CN=C1)(C2=CC=C(C)C=C2)=O',
        'n-trityl': 'N1(C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C=CN=C1',
        'n-monomethoxytrityl': 'COC(C=C1)=CC=C1C(C2=CC=CC=C2)(C3=CC=CC=C3)N4C=CN=C4',
        'n-methyltrityl': 'CC(C=C1)=CC=C1C(C2=CC=CC=C2)(C3=CC=CC=C3)N4C=CN=C4',
        'n-tert-butyloxycarbonyl': 'O=C(OC(C)(C)C)N1C=CN=C1',
        'n-2,4-dimethylpent-3-yloxycarbonyl': 'O=C(OC(C(C)C)C(C)C)N1C=CN=C1',
        'n-benzyloxymethyl': '[N+]1(COCC2=CC=CC=C2)=CNC=C1',
        'n-tert-butoxymethyl': 'CC(C)(C)OC[N+]1=CNC=C1',
    }

    his_protecting_groups_removed_by_acid_smarts = {
        'n-tosyl': '[#8]=[#16](-[#7]1:[#6]:[#6]:[#7]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1)=[#8]',
        'n-trityl': '[#7]1(-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:[#7]:[#6]:1',
        'n-monomethoxytrityl': '[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'n-methyltrityl': '[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'n-tert-butyloxycarbonyl': '[#8]=[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'n-2,4-dimethylpent-3-yloxycarbonyl': '[#8]=[#6](-[#8]-[#6](-[#6](-[#6])-[#6])-[#6](-[#6])-[#6])-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'n-benzyloxymethyl': '[#7+]1(-[#6]-[#8]-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#7H]:[#6]:[#6]:1',
        'n-tert-butoxymethyl': '[#6]-[#6](-[#6])(-[#6])-[#8]-[#6]-[#7+]1:[#6]:[#7H]:[#6]:[#6]:1',
    }

    his_protecting_groups_removed_by_acid_acronym_smiles = {
        'ntos': 'O=S(N1C=CN=C1)(C2=CC=C(C)C=C2)=O',
        'ntrt': 'N1(C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C=CN=C1',
        'nmtt': 'COC(C=C1)=CC=C1C(C2=CC=CC=C2)(C3=CC=CC=C3)N4C=CN=C4',
        'nmmt': 'CC(C=C1)=CC=C1C(C2=CC=CC=C2)(C3=CC=CC=C3)N4C=CN=C4',
        'nboc': 'O=C(OC(C)(C)C)N1C=CN=C1',
        'ndoc': 'O=C(OC(C(C)C)C(C)C)N1C=CN=C1',
        'nbom': '[N+]1(COCC2=CC=CC=C2)=CNC=C1',
        'nbum': 'CC(C)(C)OC[N+]1=CNC=C1',
    }

    his_protecting_groups_removed_by_acid_acronym_smarts = {
        'ntos': '[#8]=[#16](-[#7]1:[#6]:[#6]:[#7]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1)=[#8]',
        'ntrt': '[#7]1(-[#6](-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#6]:[#7]:[#6]:1',
        'nmtt': '[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'nmmt': '[#6]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'nboc': '[#8]=[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'ndoc': '[#8]=[#6](-[#8]-[#6](-[#6](-[#6])-[#6])-[#6](-[#6])-[#6])-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'nbom': '[#7+]1(-[#6]-[#8]-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2):[#6]:[#7H]:[#6]:[#6]:1',
        'nbum': '[#6]-[#6](-[#6])(-[#6])-[#8]-[#6]-[#7+]1:[#6]:[#7H]:[#6]:[#6]:1',
    }

    his_protecting_groups_removed_by_base_smiles = {
        'N-9-fluorenylmethoxycarbonyl': 'O=C(OC1C(C=CC=C2)=C2C3=C1C=CC=C3)N4C=CN=C4',
        '2,6-dimethoxybenzoyl': 'O=C(C1=C(OC)C=CC=C1OC)N2C=CN=C2',
    }

    his_protecting_groups_removed_by_base_smarts = {
        'N-9-fluorenylmethoxycarbonyl': '[#8]=[#6](-[#8]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'N-2,6-dimethoxybenzoyl': '[#8]=[#6](-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6]:[#6]:[#6]:1-[#8]-[#6])-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
    }

    his_protecting_groups_removed_by_base_acronym_smiles = {
        'fmoc': 'O=C(OC1C(C=CC=C2)=C2C3=C1C=CC=C3)N4C=CN=C4',
        'dmbz': 'O=C(C1=C(OC)C=CC=C1OC)N2C=CN=C2',
    }

    his_protecting_groups_removed_by_base_acronym_smarts = {
        'fmoc': '[#8]=[#6](-[#8]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
        'dmbz': '[#8]=[#6](-[#6]1:[#6](-[#8]-[#6]):[#6]:[#6]:[#6]:[#6]:1-[#8]-[#6])-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
    }

    other_his_protecting_groups_smiles = {
        'N-2,4-dinitrophenyl': 'O=C(C1=C([N+]([O-])=O)C=C([N+]([O-])=O)C=C1)N2C=CN=C2',
    }

    other_his_protecting_groups_smarts = {
        'N-2,4-dinitrophenyl': '[#8]=[#6](-[#6]1:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',

    }

    other_his_protecting_groups_acronym_smiles = {
        'dnp': 'O=C(C1=C([N+]([O-])=O)C=C([N+]([O-])=O)C=C1)N2C=CN=C2'
    }

    other_his_protecting_groups_acronym_smarts = {
        'dnp': '[#8]=[#6](-[#6]1:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1)-[#7]1:[#6]:[#6]:[#7]:[#6]:1',
    }

    ser_thr_hyp_protecting_groups_removed_by_acid_smiles = {
        'benzyl': 'CC1=CC=CC=C1',
        'cyclohexyl;': 'C1CCCCC1',
        'tert-butyl': 'CC(C)C',
        'trityl': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        'tert-butyldimethylsilyl': 'CSi(C(C)(C)C)C',
        'pseudoprolines': 'CC1(C)NC(C(O)=O)CO1'
    }

    ser_thr_hyp_protecting_groups_removed_by_acid_smarts = {
        'benzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'cyclohexyl;': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        'tert-butyl': '[#6]-[#6](-[#6])-[#6]',
        'trityl': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tert-butyldimethylsilyl': '[#6]-[Si](-[#6](-[#6])(-[#6])-[#6])-[#6]',
        'pseudoprolines': '[#6]-[#6]1(-[#6])-[#7]-[#6](-[#6](-[#8])=[#8])-[#6]-[#8]-1',
    }

    ser_thr_hyp_protecting_groups_removed_by_acid_acronym_smiles = {
        'bn': 'CC1=CC=CC=C1',
        'chx': 'C1CCCCC1',
        'bu': 'CC(C)C',
        'trt': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        'tbdms': 'C[Si](C(C)(C)C)C',
    }

    ser_thr_hyp_protecting_groups_acronym_removed_by_acid_smarts = {
        'bn': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'chx': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        'bu': '[#6]-[#6](-[#6])-[#6]',
        'trt': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tbdms': '[#6]-[Si](-[#6](-[#6])(-[#6])-[#6])-[#6]',
    }

    other_ser_thr_hyp_protecting_groups_smiles = {
        'tert-butyldiphenylsilyl': 'CC(Si(C1=CC=CC=C1)C2=CC=CC=C2)(C)C',
        '4,5-dimethoxy-2-nitrobenzyloxycarbonyl': 'CC1=CC(OC)=C(OC)C=C1[N+]([O-])=O',
        'propargyloxycarbonyl': 'C#COC(C)=O',
    }

    other_ser_thr_hyp_protecting_groups_smarts = {
        'tert-butyldiphenylsilyl': '[#6]-[#6](-[Si](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6])-[#6]',
        '4,5-dimethoxy-2-nitrobenzyloxycarbonyl': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'propargyloxycarbonyl': '[#6]#[#6]-[#8]-[#6](-[#6])=[#8]',
    }

    other_ser_thr_hyp_protecting_groups_acronym_smiles = {
        'tbdps': 'CC([Si](C1=CC=CC=C1)C2=CC=CC=C2)(C)C',
        'dmnb': 'CC1=CC(OC)=C(OC)C=C1[N+]([O-])=O',
        'poc': 'C#COC(C)=O'
    }

    other_ser_thr_hyp_protecting_groups_acronym_smarts = {
        'tbdps': '[#6]-[#6](-[Si](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6])-[#6]',
        'dmnb': '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#8]-[#6]):[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'poc': '[#6]#[#6]-[#8]-[#6](-[#6])=[#8]',
    }

    tyr_protecting_groups_removed_by_acid_smiles = {
        'benzyl': 'CC1=CC=CC=C1',
        'tert-butyl': 'CC(C)C',
        '2,6-dichlorobenzyl': 'CC1=C(Cl)C=CC=C1Cl',
        '2-bromobenzyl': 'CC1=CC=CC=C1Br',
        'benzyloxycarbonyl': 'O=COCC1=CC=CC=C1',
        '2-bromobenzyloxycarbonyl': 'O=COCC1=CC=CC=C1Br',
        '3-pentyl': 'CCCCC',
        'tert-butyloxycarbonyl': 'O=COC(C)(C)C',
        'trityl': 'C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3',
        '2-chlorotrityl': 'ClC1=CC=CC=C1C(C2=CC=CC=C2)C3=CC=CC=C3',
        'tert-butyldimethylsilyl': 'C[Si](C(C)(C)C)C',
        '4-(3,6,9-trioxadecyl)oxybenzyl': 'CC1=CC=C(OCCOCCOCCOC)C=C1'
    }

    tyr_protecting_groups_removed_by_acid_smarts = {
        'benzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tert-butyl': '[#6]-[#6](-[#6])-[#6]',
        '2,6-dichlorobenzyl': '[#6]-[#6]1:[#6](-[#17]):[#6]:[#6]:[#6]:[#6]:1-[#17]',
        '2-bromobenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#35]',
        'benzyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '2-bromobenzyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#35]',
        '3-pentyl': '[#6]-[#6]-[#6]-[#6]-[#6]',
        'tert-butyloxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        'trityl': '[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)(-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        '2-chlorotrityl': '[#17]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tert-butyldimethylsilyl': '[#6]-[Si](-[#6](-[#6])(-[#6])-[#6])-[#6]',
        '4-(3,6,9-trioxadecyl)oxybenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]):[#6]:[#6]:1',
    }

    tyr_protecting_groups_removed_by_acid_acronym_smiles = {
        'bn': 'CC1=CC=CC=C1',
        'bu': 'CC(C)C',
        'dcb': 'CC1=C(Cl)C=CC=C1Cl',
        'brbn': 'CC1=CC=CC=C1Br',
        'z': 'O=COCC1=CC=CC=C1',
        'brz': 'O=COCC1=CC=CC=C1Br',
        'pen': 'CCCCC',
        'boc': 'O=COC(C)(C)C',
        '2-cl-trt': 'ClC1=CC=CC=C1C(C2=CC=CC=C2)C3=CC=CC=C3',
        'tbdms': 'C[Si](C(C)(C)C)C',
        'tegb': 'CC1=CC=C(OCCOCCOCCOC)C=C1',
    }

    tyr_protecting_groups_removed_by_acid_acronym_smarts = {
        'bn': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'bu': '[#6]-[#6](-[#6])-[#6]',
        'dcb': '[#6]-[#6]1:[#6](-[#17]):[#6]:[#6]:[#6]:[#6]:1-[#17]',
        'brbn': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#35]',
        'z': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'brz': '[#8]=[#6]-[#8]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#35]',
        'pen': '[#6]-[#6]-[#6]-[#6]-[#6]',
        'boc': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        '2-cl-trt': '[#17]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
        'tbdms': '[#6]-[Si](-[#6](-[#6])(-[#6])-[#6])-[#6]',
        'tegb': '[#6]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]):[#6]:[#6]:1',
    }

    other_tyr_protecting_group_smiles = {
        'allyl': 'CC=C',
        'o-nitrobenzyl': 'CC1=CC=CC=C1[N+]([O-])=O',
        'propargyloxycarbonyl': 'C#COC(C)=O',
        'boc-n-methyl-n-[2-(methylamino)ethyl]carbamoyl': 'O=C(C(OC(C)(C)C)=O)N(C)CCNC'
    }

    other_tyr_protecting_group_smarts = {
        'allyl': '[#6]-[#6]=[#6]',
        'o-nitrobenzyl': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'propargyloxycarbonyl': '[#6]#[#6]-[#8]-[#6](-[#6])=[#8]',
        'boc-n-methyl-n-[2-(methylamino)ethyl]carbamoyl': '[#8]=[#6](-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])=[#8])-[#7](-[#6])-[#6]-[#6]-[#7]-[#6]',
    }

    other_tyr_protecting_group_acronym_smiles = {
        'al': 'CC=C',
        'onb': 'CC1=CC=CC=C1[N+]([O-])=O',
        'poc': 'C#COC(C)=O',
        'boc-nmec': 'O=C(C(OC(C)(C)C)=O)N(C)CCNC',
    }

    other_tyr_protecting_group_acronym_smarts = {
        'al': '[#6]-[#6]=[#6]',
        'onb': '[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]',
        'poc': '[#6]#[#6]-[#8]-[#6](-[#6])=[#8]',
        'boc-nmec': '[#8]=[#6](-[#6](-[#8]-[#6](-[#6])(-[#6])-[#6])=[#8])-[#7](-[#6])-[#6]-[#6]-[#7]-[#6]',
    }

    trp_protecting_groups_removed_by_acid_smiles = {
        'formyl': 'C=O',
        'tert-butyloxycarbonyl': 'O=COC(C)(C)C',
        'cyclohexyloxycarbonyl': 'O=COC1CCCCC1',
        'mesityl-2-sulfonyl': 'CC1=C(S(=O)(NC(N)=N)=O)C(C)=CC(C)=C1',
        'allyloxycarbonyl': 'O=COCC=C',
    }

    trp_protecting_groups_removed_by_acid_smarts = {
        'formyl': '[#6]=[#8]',
        'tert-butyloxycarbonyl': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        'cyclohexyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        'mesityl-2-sulfonyl': '[#6]-[#6]1:[#6](-[#16](=[#8])(-[#7]-[#6](-[#7])=[#7])=[#8]):[#6](-[#6]):[#6]:[#6](-[#6]):[#6]:1',
        'allyloxycarbonyl': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
    }

    trp_protecting_groups_removed_by_acid_acronym_smiles = {
        'for': 'C=O',
        'boc': 'O=COC(C)(C)C',
        'hoc': 'O=COC1CCCCC1',
        'mts': 'CC1=C(S(=O)(NC(N)=N)=O)C(C)=CC(C)=C1',
        'alloc': 'O=COCC=C'
    }

    trp_protecting_groups_removed_by_acid_acronym_smarts = {
        'for': '[#6]=[#8]',
        'boc': '[#8]=[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]',
        'hoc': '[#8]=[#6]-[#8]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
        'mts': '[#6]-[#6]1:[#6](-[#16](=[#8])(-[#7]-[#6](-[#7])=[#7])=[#8]):[#6](-[#6]):[#6]:[#6](-[#6]):[#6]:1',
        'alloc': '[#8]=[#6]-[#8]-[#6]-[#6]=[#6]',
    }

    local = list(locals().values())
    keys = list(locals().keys())

    total = 0
    for i in range(0, len(local)):
        if 'smiles' in keys[i]:
            total += len(local[i])
            print (keys[i])

if __name__ == '__main__':
    _get_common_amino_acid_protecting_groups()