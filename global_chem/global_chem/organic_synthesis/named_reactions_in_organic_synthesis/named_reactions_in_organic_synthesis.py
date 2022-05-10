#!/usr/bin/env python3
#
# GlobalChem - Named Reactions in Organic Synthesis
#
# -------------------------------------------------

class NamedReactionsInOrganicSynthesis(object):
  
  def __init__(self):
    
      self.name = 'named_reactions_in_organic_synthesis'
      
  @staticmethod
  def get_smiles():

      smiles = {
        'lieben_haloform_reaction': 'CC(C)=O.CC([O-])=O',
        'benzilic_acid_rearrangement': 'O=C(C1=CC=CC=C1)C(C2=CC=CC=C2)=O.OC(C(O)=O)(C3=CC=CC=C3)C4=CC=CC=C4',
        'aldol_reaction': 'CC(C(C)=O)C.CC(C(C)=O)(C(O)C)C',
        'dieckmann_condensation': 'O=C(O)CCCCC(O)=O.O=C1C(C(O)=O)CCC1',
        'strecker_reation': 'CC=O.O=C(O)CN',
        'hofmann_elimination': 'CCC(C)N.CCC=C',
        'wiallmson_ether_synthesis': 'CO.C[F,Cl,Br,I].COC',
        'cannizzaro_reaction': 'CC=O.CCO',	# R2 CC(O)=O
        'wurtz_coulping': 'CBr.CC',
        'kolbe-schmitt_reaction': 'OC2=CC=CC=C2.OC1=C(C(O)=O)C=CC=C1',
        'pinacol_rearrangement': 'OC(C)(C)C(O)(C)C.O=C(C)C(C)(C)C',
        'semipinacol_rearrangement': 'OC(C)(C)C([F,Cl,Br,I])(C)C.O=C(C)C(C)(C)C',
        'acyloin_condensation': 'COC(C)=O.CC(O)C(C)=O',
        'hunsdiecker_reaction_borodin_reaction': 'CC([O-])=O.CBr',
        'perkin_reaction': r'O=CC1=CC=CC=C1.O=C(C)OC(C)=O.O=C(O)/C=C\C2=CC=CC=C2',
        'glaser_coupling_reaction': 'CC#C[H].CC#CC#CC',
        'lossen_rearrangement': 'CC(N([H])OC(C)=O)=O.CN=C=O',
        'reimer-tiemann_reaction': 'OC1=CC=CC=C1.OC2=C(C([H])=O)C=CC=C2',
        'friedel-crafts_acylation': 'C1=CC=CC=C1.O=C(C)C2=CC=CC=C2',
        'friedel-crafts_alkylation': 'C1=CC=CC=C1.CC2=CC=CC=C2',
        'malonic_ester_synthesis': 'O=C(CC(OC)=O)OC.O=C(OC)CC',
        'pinner_reaction': 'CC#N.CC(OC)=N',
        'koenigs-knorr_glycosidation': 'CC(OCC1C(OC(C)=O)C(OC(C)=O)C(OC(C)=O)C(Br)O1)=O.CC(OCC2C(OC(C)=O)C(OC(C)=O)C(OC(C)=O)C(OC)O2)=O',
        'skraup_and_doebner-miller_reaction': r'NC1=CC=CC=C1.O=C/C=C/C.CC2=NC3=CC=CC=C3C=C2',
        'ciamician-dennstedt_rearrangement': 'C1=CNC=C1.[Cl,Br,I]C2=CC=CN=C2',
        'fries-_photo-fries_rearrangement': 'O=C(C)OC1=CC=CC=C1.OC2=CC=C(C(C)=O)C=C2',
	'anionic_ortho-fries_rearrangement': 'O=C(C)OC1=CC=CC=C1.O=C(C)C2=C(O)C=CC=C2',
        'hell-volhard-zelinsky_reaction': 'CCC(O)=O.CC([Cl,Br])C(O)=O',
        'hofmann_rearrangement': 'CC(N)=O.CN',
        'hantzsch_dihydropyridine_synthesis': 'CC(CC(C)=O)=O.CC=O.CC1=C(C(C)=O)C(C)C(C(C)=O)=C(C)N1',
        'combes_quinoline_synthesis': 'NC1=CC=CC=C1.CC(CC(C)=O)=O.CC2=CC(C)=NC3=CC=CC=C32',
        'fischer_indole_synthesis': 'NNC1=CC=CC=C1.CCC(C)=O.CC2=C(C)NC3=CC=CC=C32',
        'hofmann–löffler-freytag_reaction': 'CCCCCN([Cl,Br,I])C.CC1CCCN1C',
        'michael_addition': 'CC(CC(OCC)=O)=O.C=CC(C)=O.CC(C(C(OCC)=O)CCC(C)=O)=O',
        'von_pechmann_reaction': 'OC1=CC=CC=C1.CC(CC(OCC)=O)=O.CC(C2=CC=CC=C2O3)=CC3=O',
        'paal-knorr_furan_synthesis': 'CC(CCC(C)=O)=O.CC1=CC=C(C)O1',
        'paal-knorr_pyrrole_synthesis': 'CC(CCC(C)=O)=O.CN.CC1=CC=C(C)N1C',
        'sandmeyer_reaction': 'NC1=CC=CC=C1.[F,Cl,Br,I]C2=CC=CC=C2',
        'schotten-baumann_reaction': 'CC(Cl)=O.CN.CC(NC)=O',	# from amine, R CC(OC)=O from alcohol 
        'buchner_method_of_ring_enlargement_buchner_reaction': 'C1=CC=CC=C1.O=C(OCC)C=[N+]=[N-].O=C(OCC)C2C=CC=CC=C2',
        'curtius_rearrangement': 'CC(N=[N+]=[N-])=O.CN',	# from water, R CNC(=O)OC from alcohol, R CNC(=O)NC from amine 
        'beckman_rearrangement': 'CC(=NO)C.CC(NC)=O',
        'knorr_pyrrole_synthesis': 'O=C(C(C)N)C.O=C(OCC)CC(C)=O.CC1=C(C(OCC)=O)C(C)=C(C)N1',
        'claisen_condensation': 'CCOC(C)=O.CCOC(CC(C)=O)=O',
        'gabriel_synthesis': 'C[F,Cl,Br,I].CN',
        'japp-klingemann_reaction': r'[N+]=NC1=CC=CC=C1.CC(C(C)C(O)=O)=O.C/C(C(OC)=O)=N\NC1=CC=CC=C1',
        'reformatsky_reaction': '[F,Cl,Br,I]CC(OC)=O.C=O.COC(CCO)=O',
        'tishchenko_reaction': 'CC=O.CC(OCC)=O',
        'dimroth_rearrangement': 'NC1=C(C)N=NN1C.CC2=C(NC)NN=N2',
        'biginelli_reaction': 'NC(N)=O.O=CC1=CC=CC=C1.CC(CC(OCC)=O)=O.O=C(OCC)C(C(N2)C3=CC=CC=C3)=C(C)NC2=O',
        'darzens_glycidic_ester_condensation': 'CC(C)=O.CC(C)C(C)=O',
        'bischler-napieralski_isoquinoline_synthesis': 'CC(NCCC1=CNC2=CC=CC=C21)=O.CC3=NCCC4=C3NC5=CC=CC=C54',
        'dienone-phenol_rearrangement': 'CC(C=C1)(C)C=CC1=O.CC2=C(C)C=C(O)C=C2',
        'pomeranz-fritsch_reaction': 'O=CC1=CC=CC=C1.NCC(OC)OC.C23=CC=CC=C2C=NC=C3',
        'stobbe_condensation': r'CCOC(CCC(OCC)=O)=O.CC(C)=O.C/C(C)=C(CC(O)=O)\C(OCC)=O',
        'favorskii_rearrangement_and_homo-favorskii_rearrangement': '[F,Cl,Br,I]C1CCCCC1=O.O=C(O)C2CCCC2',
        'knoevenagel_condensation': r'O=CC1=CC=CC=C1.COC(CC(OC)=O)=O.O=C(OC)/C(C(OC)=O)=C/C2=CC=CC=C2',
        'nef_reaction': 'O=[N]([O])C(C)C.O=C(C)C',
        'smiles_rearrangement': '[F,Cl,Br,I]C1=C(OCCN)C=CC([N+]([O-])=O)=C1.[F,Cl,Br,I]C2=C(NCCO)C=CC([N+]([O-])=O)=C2',	# X-Y rearrangment where X(here O) and Y(here NH2) can be many different groups, also ortho- and para- substituents(here halogens and NO2) can be different groups
        'wacker_oxidation': 'C=C.CC=O',
        'henry_reaction': 'CC[N+]([O-])=O.CC(C)=O.C/C([N+]([O-])=O)=C(C)\C',	# from dehydration, CC(C(C)=C)[N+]([O-])=O from oxidation, CC(C(O)(C)C)N from reduction
        'michaelis-arbuzov_reaction': 'CP(C)OC.CP(C)(C)=O',	# phosphinite to phosphine, CP(OC)OC.CP(C)(OC)=O phosphonite to phosphinate, COP(OC)OC.CP(OC)(OC)=O phosphite ester to phosphonate oxide
        'gattermann_and_gattermann-koch_formylation': 'C#N.O=CC1=CC=CC=C1',
        'chugaev_elimination_xanthate_ester_pyrolysis': r'CC(C)C(C)(C)O.C/C(C)=C(C)/C',
        'baeyer-villiger_oxidation': 'CC(C)=O.CC(OO)=O.CC(OC)=O',
        'barbier_coupling_reaction': 'C[Cl,Br,I].CC(C)=O.CC(O)(C)C',
        'prins_reaction': 'CC=C.CC(C)=O.OC(C)CCO',	# to 1,3 diol, C/C=C/CO allylic alcohol in absence of water, CC1OCOCC1 dioxane with excess of CC(C)=O
        'wagner-meerwein_rearrangement': r'CCC[F,Cl,Br,I].C/C=C/C',	# or to CC([F,Cl,Br,I])C (hydrid shift)
        'grignard_reaction': 'C=O.C[Mg]Br.CCO',	# to primary alcohol, other aldehyds to secondary alcohols CC=O.CC(O)C, ketones to tertiary alcohols CC(C)=O.CC(C)(O)C,…  
        'demjanov_rearrangement': 'NCC1CCC1.OC2CCCC2',	# ring expansion, or ring reduction NC1CCCC1.OCC2CCC2
	'tiffeneau-demjanov_rearrangement': 'OC1(CN)CCCC1.O=C2CCCCC2',
        'ullmann_reaction_coupling_biaryl_synthesis': '[F,Cl,Br,I]C1=CC=CC=C1.C2(C3=CC=CC=C3)=CC=CC=C2',
        'feist–bénary_furan_synthesis': 'CC(CC(C)=O)=O.O=C(C)C([Cl,Br,I])C.CC1=C(C(C)=O)C(C)=C(C)O1',
        'wolff_rearrangement': 'O=C(C)C(C)=[N+]=[N-].OC(C(C)C)=O',
        'benzoin_and_retro-benzoin_condensation': 'O=CC1=CC=CC=C1.O=C(C(O)C2=CC=CC=C2)C3=CC=CC=C3',
        'mannich_reaction': 'CNC.C=O.CC(C(C)=O)C.CN(CC(C)(C(C)=O)C)C',
        'nazarov_cyclization': r'C/C=C/C(/C=C/C)=O.O=C1C=C(C)C(C)C1',
        'ullmann_biaryl_ether_and_biaryl_amine_synthesis_condensation': '[F,Cl,Br,I]C1=CC=CC=C1.C2(C3=CC=CC=C3)=CC=CC=C2',
        'eschweiler-clarke_methylation': 'CNC.O=C.OC=O.CN(C)C',
        'staudinger_ketene_cycloaddition': r'C/C(C)=N\C.CC(C)=C=O.CN1C(C)(C)C(C)(C)C1=O',
        'acetoacetic_ester_synthesis': 'CC(CC(OC)=O)=O.CC(CC)=O',
        'dakin_oxidation': 'OC1=CC=C(C=O)C=C1.OC2=CC=C(C(C)=O)C=C2.OC3=CC=C(O)C=C3',
        'paterno-büchi_reaction': r'CC(C)=O.C/C(C)=C(C)\C.CC1(C)C(C)(C)C(C)(C)O1',	# stereochemistry when R groups are different on the 2nd reactant
        'prilezhaev_reaction': r'C/C(C)=C(C)\C.CC(OO)=O.CC1(C)OC1(C)C',
        'pummerer_rearrangement': 'CS(CC)=O.CC(OC(C)=O)=O.CSC(OC(C)=O)C',
        'finkelstein_reaction': 'C[F,Br,I].[Cl-].CCl',	# halogen replacement with any combination 
        'regitz_diazo-transfer_reaction': 'CC(CC)=O.[N]=[N+]=[N-].CC(C(C)=[N+]=[N-])=O',
        'pictet-spengler_tetrahydroisoquinoline_synthesis': 'NCCC1=CC=CC=C1.CC=O.CC2C3=CC=CC=C3CCN2',	# from Phenylethylamine, NCCC1=CNC2=CC=CC=C21.CC=O.CC(NCC3)C4=C3C5=CC=CC=C5N4 from Tryptamine
        'wolff-kishner_reduction': 'CC(C)=O.CCC',
        'madelung_indole_synthesis': 'CCC1=CC=CC=C1NC(C)=O.CC2=C(C)NC3=CC=CC=C32',
        'claisen_rearrangement': r'C=CCO/C(C)=C/C.O=C(C(CC=C)C)C',
        'clemmensen_reduction': 'CC(C)=O.CCC',	# reduction with Zn		
        'wharton_olefin_synthesis': r'O=C(C)C1OC1C.C/C=C\C(C)O',
        'chichibabin_amination_reaction': 'CC1=CC=CC=N1.CC2=CC=CC(N)=N2',
        'ferrier_rearrangement_reaction': 'O=C(C)OC1C(COC(C)=O)OC=CC1OC(C)=O.O=C(C)OC2C(COC(C)=O)OC(OC)C=C2',	
	'ferrier_II_carbocyclization_reaction': 'C=C1C(OC)C(OC)C(OC)C(OC)O1.O=C2C(OC)C(OC)C(OC)C(O)C2', 
        'houben-hoesch_reaction': 'CC#N.OC1=CC(O)=CC=C1.OC2=CC(O)=CC=C2C(C)=O',
        'aza-wittig_reaction': r'CP(C)(C)=NC.CC(C)=O.C/C(C)=N/C',	# with ketone to imine, with CO2 to isocyanate, with CS2 to isothiocyanate, with H2O to primary amine (=staudinger),…
        'meisenheimer_rearrangement': '[O-][N+](C)(C)C.CN(OC)C',	# [1,2] rearrangement, [O-][N+](C)(C(C)C=C)C.CN(OC/C=C/C)[2,3] rearrangement
        'staudinger_reaction': 'CP(C)(C)=NC.CN',	# = aza-wittig with H2O
        'wittig_reaction': r'C[P+](C)([C-](C)C)C.CC(C)=O.C/C(C)=C(C)\C',
        'wohl-ziegler_bromination': 'C1CCC=CC1.BrC2C=CCCC2',
        'passerini_multicomponent_reaction': 'OC(C)=O.CC(C)=O.[C-]#[N+]C.CC(OC(C)(C)C(N(C)C)=O)=O',
        'meyer-schuster_and_rupe_rearrangement': r'CC(CC)(C#C)O.C/C(CC)=C\C=O',	# mayer-schuster (1,3 shift of hydroxyl group), rupe (1,2 shift of hydroxyl group) to C/C(C(C)=O)=C/C 
        'schmidt_reaction': 'CC(C)=O.[H][N-][N+]#N.CC(N(C))=O',	# with ketones to amides, with carboxylic acids to amines (CC(O)=O.[H][N-][N+]#N.CN)
        'amadori_reaction_rearrangement': r'OC(/C=N/C)C.CC(CNC)=O',
        'meerwein-ponndorf-verley_reduction': 'CC(C)=O.OC(C)C',	# reverse reaction of oppenauer_oxidation
        'stephen_aldehyde_synthesis': 'CC#N.CC=O',
        'diels-alder_cycloaddition': 'C=CC=C.C=C.C1=CCCCC1',
        'neber_rearrangement': 'CC(CC)=O.CC(C(N)C)=O',
        'balz-schiemann_reaction': 'NC1=CC=CC=C1.FC2=CC=CC=C2',
        'dakin-west_reaction': 'OC(C(C)N)=O.CC(C(C)NC(C)=O)=O',
        'stevens_rearrangement': 'CC[N+](C)(C)C.CC(N(C)C)C',	# quaternary ammonium salts to amines or CC[S+](C)C.CC(SC)C sulfonium salts to sulfides 
        'nenitzescu_indole_synthesis': r'O=C1C=CC(C=C1)=O.O=C(OC)/C=C(C)/NC.OC2=CC=C(N(C)C(C)=C3C(OC)=O)C3=C2',
        'criegee_oxidation': 'OCCO.C=O',
        'riley_selenium_dioxide_oxidation': 'CCC(CC)=O.CCC(C(C)=O)=O',	# R2: CC(C(CC)=O)=O oxidation of methylene groups adjacent to carbonyls
        'baker-venkataraman_rearrangement': 'CC(OC1=CC=CC=C1C(CC)=O)=O.OC2=CC=CC=C2C(C(C)C(C)=O)=O',
        'prévost_reaction': r'C/C=C\C.OC(C(C)O)C',
        'arndt-eistert_homologation': 'CC(O)=O.OC(CC)=O',
        'payne_rearrangement': 'CC1C(CO)O1.CC(O)C2OC2',
        'robinson_annulation': r'O=C(C(C)C)C(C)C.O=C(/C(C)=C(C)/C)CC.O=C1C(C)C(C)(C)C(C)(C)C(C(C)C)=C1C',
        'claisen-ireland_rearrangement': 'C=CCOC(C)=O.OC(CCC=C)=O',
        'oppenauer_oxidation': 'OC(C)C.CC(C)=O', # reverse reaction of meerwein-ponndorf-verley_reduction
        'sommelet-hauser_rearrangement': 'C[N+](C)(C)CC1=CC=CC=C1.CC2=CC=CC=C2CN(C)C',
        'meerwein_arylation': r'N#[N+]C1=CC=CC=C1.C/C(C2=CC=CC=C2)=C(C)\[N+]([O-])=O.C/C=C([N+]([O-])=O)/C', # instead of NO2 also other EWGs possible
        'quasi-favorskii_rearrangement': 'CC([Cl,Br])C(C(C)(C)C)=O.CC(C(C)(C)C)C(OC)=O',
        'snieckus_directed_ortho_metalation': 'O=C(N(C)C)C1=CC=CC=C1.O=C(N(C)C)C2=C([Li])C=CC=C2',
        'carroll_rerrangement': 'CC(CC(OCC=C)=O)=O.CC(CCC=C)=O',
        'cope_rearrangement': r'C=CCC(C)C=C.C=CCC/C=C/C',
        'ramberg–bäcklund_rearrangement': r'CCS(C(Br)C)(=O)=O.C/C=C\C',
        'wittig-[1,2]-_and_[2,3]-rearrangement': r'CCOC(CC)/C=C/C.CC(C(/C=C/C)CC)O', # [2,3] product CC(C(/C=C/CC)C)O
        'alder_reaction': '', # alder without diels? haha
        'hetero_diels-alder_reaction': 'C=CC=C.C1=CCCCC1',
        'birch_reduction': 'C1=CC=CC=C1.C2=CCC=CC2',
        'jones_oxidation_primary_alcohol': 'CCO.CC(O)=O',
	'jones_oxidation_secondary_alcohol': 'CC(O)C.CC(C)=O',
        'peterson_olefination': r'CC(CC)=O.C[C-]([H])[Si](C)(C)C.C/C(CC)=C\C',	# R2 C/C(CC)=C/C (cis/trans alkenes)
        'ritter_reaction': 'CC(C)(O)C.N#CC.CC(C)(NC(C)=O)C',
        'cope_elimination': 'CCCN(C)C.CN(O)C',	# R2 C=CC
        'cornforth_rearrangement': 'CC1=NC(C(C)=O)=C(CC)O1.CC2=NC(C(CC)=O)=C(C)O2',
        'bamford-stevens-shapiro_olefination': r'CC(CC)=O.C/C=C/C',	# R2 C/C=C\C (cis/trans alkenes)
        'grob_fragmentation': 'BrC1CCC(Br)CC1.C=CCCC=C',	# original reaction or [F,Cl,Br,I]CCCOC.C[O+]=C + R2 C=C (similar to wharton fragmentation)
        'wharton_fragmentation': 'OC(C)(C)C(C)(C)C(C)(C)[F,Cl,Br,I].O=C(C)C',	# R2 C/C(C)=C(C)\C
        'stork_enamine_synthesis': 'O=C1CCCCC1.C2CNCC2.C3(N4CCCC4)=CCCCC3',	# enamine = temporary activator, further reactions: alkylation, acylation, and conjugate addition	
        'alkene_olefin_metathesis': r'CC/C=C\CC.C/C=C\C.C/C=C\CC',	# name reaction? haha 
        'brown_hydroboration_reaction': 'CC=C.CCCO',	# with BH3
        'kornblum_oxidation': '[Br,I]CC1=CC=CC=C1.C=CC2=CC=CC=C2',	# or with carbonyl group [Br,I]CC1=CC=CC=C1.C=CC2=CC=CC=C2
        'brook_rearrangement': 'C[Si](C)(CO)C.C[Si](C)(OC)C',
        'doering-laflamme_allene_synthesis': r'C/C(C)=C(C)/C.[Cl,Br]C[Cl,Br].CC(C)=C=C(C)C',
        'horner-wadsworth-emmons_olefination': r'O=P(OC)(CC)OC.CC=O.C/C=C/C',
        'simmons-smith_cyclopropanation': r'C/C(C)=C(C)/C.ICI.CC1(C)C(C)(C)C1',
        'heine_reaction': 'O=C(C)N1CC1.CC2=NCCO2',  # Aziridine synthesis into oxazolines are important, transformations for potential warheads this one makes me happy.
        'ugi_multicomponent_reaction': 'CC(O)=O.CN.CC(C)=O.[C-]#[N+]C.CC(N(C)C(C)(C)C(NC)=O)=O',
        'vilsmeier-haack_formylation': 'CN(C)C1=CC=CC=C1.O=CN(C)C.O=CC2=CC=C(N(C)C)C=C2',
        'vinylcyclopropane-cyclopetene_rearrangement': 'C=CC1CC1.C2=CCCC2',
        'barton_nitrile_ester_reaction': 'CCCCON=O.OCCCCN=O',
        'kröhnke_pyridine_synthesis': 'O=CC[N+]1=CC=CC=C1.C2=CN=CC=C2', # hahaha
        'barton_radical_decarboxylation_reaction': 'CC(O)=O.C',	# (R2 O=C=O) or with acyl chloride CC(Cl)=O.C
        'corey-chaykovsky_epoxidation_and_cyclopropanation': 'CC(C)=S.CC(C)=O.CC1(C)OC1',
        'demayo_cycloaddition_enone-alkene_[2+2]_photocycloaddition': 'CC(CC(C)=O)=O.C/C=C\C.CC(C(C)C(C)CC(C)=O)=O',	# 1,3 diketone with alkene to cyclobutane (= [2+2] cycloaddition) further to 1,5 diketone   
        'nagata_hydrocyanation_reaction': r'CC(/C=C/C)=O.CC(CC(C#N)C)=O',
        'castro-stevens_coupling': '[Cl,Br,I]C1=CC=CC=C1.[Cu]C#CC.CC#CC2=CC=CC=C2',
        'corey-winter_olefination': r'CC(C)(O)C(C)(C)O.C/C(C)=C(C)/C',
        'pfitzner-moffatt_oxidation': 'CCO.CC=O',	# primary alcohol to aldehy (no further oxidation to carboxylic acids); secondary alcohols to ketones CC(C)O.CC(C)=O
        'eschenmoser-claisen_rearrangement': 'CCOC=C.C=CCO.O=CCCC=C',
        'oxy-cope_rearrangement': 'C=CCC(O)C=C.C=CCCCC=O',
        'tsuji-trost_reaction': r'C/C=C/CCl.C/C=C/C/C(O)=C/C',	# instead of Cl other leaving groups (carbonates/phenols/phosphates/…), instead of enolate other nucleophiles (malonates/amines/sulfones/…)
        'tsuji-wilkinson_decarbonylation': '',
        'wittig_reaction-schlosser_modification': '',
        'aza-claisen_rearrangement': '',
        'aza-cope_rearrangement': '',
        'eschenmoser-tanabe_fragmentation': '',
        'krapcho_dealkoxycarbonylation': '',
        'mitsunobo_reaction': '',
        'seyferth-gilber_homologation': '',
        'baylis-hillman_reaction': '',
        'heck_reaction': r'CC=C.C/C=C/C',
        'minisci_reaction': '',
        'mislow-evans_rearrangement': '',
        'prins-pinacol_rearrangement': '',
        'schwartz_hydrozirconation': '',
        'burgess_dehydration_reaction': '',
        'johnson-claisen_rearrangement': '',
        'aza-[2,3]-wittig_rearrangement': '',
        'corey-kim_oxidation': '',
        'eschenmoser_methenylation': '',
        'hajos-parrish_reaction': '',
        'nicholas_reaction': '',
        'bergman_cycloaromatization_reaction': '',
        'corey-fuch_alkyne_synthesis': '',
        'kumada_cross_coupling_reaction': '',
        'mcmurry_coupling': '',
        'saegusa_oxidation': '',
        'julia-lythgoe_olefination': '',
        'mukaiyama_aldol_reaction': '',
        'pauson-khand_reaction': '',
        'pinnick_oxidation': '',
        'polonovski_reaction': '',
        'stetter_reaction': '',
        'overman_rearrangement': '',
        'alkyne_metathesis': '', # Huh....
        'corey-nicolau_macrolactonization': '',
        'danishefskys_diene_cycloaddition': '',
        'rubottom_oxidation': '',
        'swern_oxidation': 'CC(C)O.CC(C)=O',
        'barton-mccombie_radical_deoxygenation_reaction': '',
        'dötz_benzannulation_reaction': '',
        'sonogashira_cross-coupling': '',
        'enders_samp_ramp_hydrazone_alkylation': '',
        'negishi_cross-coupling': '',
        'sakurai_allylation': '',
        'stille_cross-coupling': '',
        'tebbe_olefination': '',
        'davis_oxaziridine_oxidation': '',
        'nozaki-hiyama-kishi_coupling': '',
        'bartoli_indole_synthesis': '',
        'luche_reduction': '',
        'roush_asymmetric_allylation': '',
        'midland_alpine_borane_reduction': '',
        'suzuki_cross-coupling': '',
        'yamaguchi_macrolactonization': '',
        'kagan-molander_samarium-diiodide_coupling': '',
        'noyori_asymmetric_hydrogenation': '',
        'sharpless_asymmetric_dihydroxylation_reaction': '',
        'Sharpless_asymmetric_epoxidation_reaction': '',
        'corey-bakashi-shibata_reduction': '',
        'danheiser_cyclopentene_annulation': '',
        'evans_aldol_reaction': '',
        'weinreb_ketone_synthesis': '',
        'buchwald-hartwig_cross-coupling': '',
        'dess-martin_oxidation': '',
        'fleming-tamao_oxidation': '',
        'danheiser_benzannulation': '',
        'stille_carbonylative_cross-couping': '',
        'enyne_metathesis': '',
        'keck_macrolactonization': '',
        'ley_oxidation': '',
        'takai-utimoto_olefination': '',
        'stille-kelly_coupling': '',
        'kahne_glycosidation': '',
        'kulinkovich_reaction': '',
        'jacobsen-katsuki_epoxidation': '',
        'larock_indole_synthesis': '',
        'keck_asymmetric_allylation': '',
        'keck_radical_allylation': '',
        'petasis_boronic_acid-mannich_reaction': '',
        'myers_asymmetric_alkylation': '',
        'smith-tietze_multicomponent_dithiane_linchpin_coupling': '',
        'jacobsen_hydrolytic_kinetic_resolution_of_epoxides': '',
        'miyaura boration_reaction': '',
        'petasis-ferrier_rearrangement': '',
        'sharpless_asymmetric_aminohydroxylation_reaction': '',
        'shi_asymmetric_epoxidation': '',
      }
      
      return smiles

  @staticmethod
  def get_smarts():

    smarts = {

    }

    return smarts
