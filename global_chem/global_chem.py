#!/usr/bin/env python3
#
# GlobalChem - Content Variable Store
#
# -----------------------------------


class GlobalChem(object):

    __version__ = "0.1.0"
    __allow_update__ = True

    """

    GlobalChem will be the master class of all variables, as the content store grows we can use this as the parent class.

    """
    
    def _get_common_regex_patterns(self):
                
        regex_patterns = {
            'mol2': '^@<\w+?>\w+?\n[COMPOUND_ID]\n(.|\n)*?@<TRIPOS>SUBSTRUCTURE\n.*?\n'
        }
        
        return regex_patterns

    def _get_amino_acids(self):

        amino_acid_side_chains = {
            "alanine": "C",  "arginine": "CCCCNC(N)=N", "asparagine": "CCC(N)=O", "aspartic acid": "CC(O)=O",
            "cysteine": "CS", "glutamic acid": "CCC(O)=O", "glutamine": "CCC(N)=O", "glycine": "[H]",
            "histidine": "CC1=CNC=N1", "isoleucine": "C(CC)([H])C", "leucine": "CC(C)C", "lysine": "CCCCN",
            "methionine": "CCSC", "phenylalanine": "CC1=CC=CC=C1", "proline": "C2CCCN2", "serine": "CO",
            "threonine": "C(C)([H])O", "tryptophan": "CCC1=CNC2=C1C=CC=C2", "tyrosine": "CC1=CC=C(O)C=C1",
            "valine": "C(C)C"
        }

        return amino_acid_side_chains

    def _get_functional_groups_smiles(self):

        functional_groups_smies = {
            "1,1,1-trifluoroethane": "CC(F)(F)F",
            "1,1'-biphenyl": "C1(C2=CC=CC=C2)=CC=CC=C1",
            "1H-indene": "C1(CC=C2)=C2C=CC=C1",
            "1H-pyrrole": "[NH]1CCCC1",
            "2-butyne": "CC#CC",
            "2-ethyl-1-butanol": "CCC(CC)CO",
            "2-methylpenta-2,3-diene": "CC=C=C(C)C",
            "(E)-1,2-dimethyldiazene": "C/N=N/C",
            "N,N-dimethylacetamide": "CC(N(C)C)=O",
            "N-methylpropan-2-imine": "C/C(C)=N/C",
            "(Z)-N,N,N'-trimethylacetimidamide": "C/C(N(C)C)=N/C",
            "acetic anydride": "CC(=O)OC(=O)C",
            "acyl bromide": "C(=O)Br",
            "acyl chloride": "C(=O)Cl",
            "acyl fluoride": "C(=O)F",
            "acyl iodide": "C(=O)I",
            "aldehyde": "CC=O",
            "amide": "C(=O)N",
            "amino": "*N",
            "anthracene": "C12=CC=CC=C1C=C3C(C=CC=C3)=C2",
            "azide": "C([N-][N+]#N)",
            "benzene": "C1=CC=CC=C1",
            "benzene thiol": "C1=CC=C(C=C1)S",
            "bicyclohexyl": "C1CCCCC1C1CCCCC1",
            "bromine": "Br",
            "but-1-ene": "CCC=C",
            "but-1-yne": "CCC#C",
            "carbon dioxide": "O=C=O",
            "carboxylic acid": "C(=O)O",
            "chlorine": "Cl",
            "chloromethyl methyl ether": "COCCl",
            "cyclobutadiene": "C1=CC=C1",
            "cyclobutane": "C1CCC1",
            "cycloheptane": "C1CCCCCC1",
            "cyclohexane": "C1CCCCC1",
            "cyclohexa-1,3-diene": "C1=CCCC=C1",
            "cyclohexa-1,4-diene": "C1=CCC=CC1",
            "cyclohexene": "C=1CCCCC=1",
            "cyclopentane": "C1CCCC1",
            "cyclopenta-1,3-diene": "C1=CCC=C1",
            "cyclopropane": "C1CC1",
            "cyclopropene": "C1=CC1",
            "deuteroethane": "[2H][CH2]C",
            "dimethyl ether": "COC",
            "diethyl ether": "CCOCC",
            "diisopropyl ether": "CC(C)OC(C)C",
            "diamond": "C&1&1&1&1",
            "diazomethane": "C=[N+]=[N-]",
            "diammonium thiosulfate": "[NH4+].[NH4+].[O-]S(=O)(=O)[S-]",
            "enamine": "N",
            "ethane": "CC",
            "ethanethiol": "CCS",
            "ethanol": "CCO",
            "ethene": "C=C",
            "ether": "COC",
            "ester": "C(=O)OC",
            "fluorine": "F",
            "formaldehyde": "C=O",
            "furan": "C1OC=CC=1",
            "graphite": "C&1&1&1",
            "hydrogen cyanide": "C#N",
            "hydroxide": "[OH-]",
            "hydroxyl amine": "NO",
            "indane": "C1=CC=CC(CCC2)=C12",
            "ketone": "CC(=O)C",
            "methane": "C",
            "methanethiol": "CS",
            "methyl acetate": "CC(OC)=O",
            "methyl pyrrole": "CN1CCCC1",
            "methyl tert-butyl ether": "CC(C)(C)OC",
            "naphthalene": "C12=CC=CC=C1C=CC=C2",
            "nitro": "[N+](=O)[O-]",
            "nitromethane": "C[N+]([O-])=O",
            "pentalene": "C12=CC=CC1=CC=C2",
            "perhydroisoquinoline": "N1CC2CCCC2CC1",
            "phenol": "OC1CCCCC1",
            "phenyl": "C=1(C=CC=CC1)",
            "polystyrene": "c1ccccc1C&1&1",
            "primary alcohol": "O",
            "primary amine": "N",
            "propan-2-one": "CC(C)=O",
            "propanol": "CCC=O",
            "prop-1-ene": "CC=C",
            "prop-1-yne": "CC#C",
            "pyridine": "N1CCCCC1",
            "pyridine-n-oxide": "O=N1CCCCC1",
            "secondary amine": "NC",
            "spiro[5.5]undecane": "C12(CCCCC1)CCCCC2",
            "sulfoxide": "S(=O)(=O)",
            "tetramethylammonium": "C[N+](C)(C)C",
            "thiol": "S",
            "thiosulfate": "OS(=O)(=S)O",
            "trimethylamine": "CN(C)C",
            "triphenylene": "C1(C=CC=C2)=C2C(C=CC=C3)=C3C4=C1C=CC=C4",
        }

        return functional_groups_smies

    def _get_functional_groups_smarts(self):

        functional_groups_smarts = {
            "adamantane": "[#6]12-[#6]-[#6]3-[#6]-[#6](-[#6]-1)-[#6]-[#6](-[#6]-3)-[#6]-2",
            "acetic anydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
            "acetylenic carbon": "[$([CX2]#C)]",
            "acyl bromide": "[CX3](=[OX1])[Br]",
            "acyl chloride": "[CX3](=[OX1])[Cl]",
            "acyl fluoride": "[CX3](=[OX1])[F]",
            "acyl iodide": "[CX3](=[OX1])[I]",
            "aldehyde": "[CX3H1](=O)[#6]",
            "alkane": "[CX4]",
            "allenic carbon": "[$([CX2](=C)=C)]",
            "amide": "[NX3][CX3](=[OX1])[#6]",
            "amidium": "[NX3][CX3]=[NX3+]",
            "amino acid": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]",
            "azepane": "[#6]1-[#7]-[#6]-[#6]-[#6]-[#6]-[#6]-1",
            "azide": "[$(-[NX2-]-[NX2+]#[NX1]),$(-[NX2]=[NX2+]=[NX1-])]",
            "azo nitrogen": "[NX2]=N",
            "azole": "[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]",
            "azoxy nitrogen": "[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]",
            "diazene": "[NX2]=[NX2]",
            "diazo nitrogen": "[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]",
            "benzofuran": "[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#8]:2",
            "bromine": "[Br]",
            "butane": "[#6]-[#6]-[#6]-[#6]",
            "carbamate": "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]",
            "carbamic ester": "[NX3][CX3](=[OX1])[OX2H0]",
            "carbamic acid": "[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]",
            "carbo azosulfone": "[SX4](C)(C)(=O)=N",
            "carbo thiocarboxylate": "[S-][CX3](=S)[#6]",
            "carbo thioester": "S([#6])[CX3](=O)[#6]",
            "carboxylate ion": "[CX3](=O)[O-]",
            "carbonic acid": "[CX3](=[OX1])(O)O",
            "carbonic ester": "C[OX2][CX3](=[OX1])[OX2]C",
            "carbonyl group": "[$([CX3]=[OX1]),$([CX3+]-[OX1-])]",
            "carbonyl with carbon": "[CX3](=[OX1])C",
            "carbonyl with nitrogen": "[OX1]=CN",
            "carbonyl with oxygen": "[CX3](=[OX1])O",
            "carboxylic acid": "[CX3](=O)[OX1H0-,OX2H1]",
            "chlorine": "[Cl]",
            "cyanamide": "[NX3][CX2]#[NX1]",
            "cyclohexene": "[#6]1=[#6]-[#6]-[#6]-[#6]-[#6]-1",
            "cyclopentanone": "[#8]=[#6]1-[#6]-[#6]-[#6]-[#6]-1",
            "dioxolane": "[#8]1-[#6]-[#6]-[#8]-[#6]-1",
            "dithiolane": "[#16]1-[#6]-[#6]-[#16]-[#6]-1",
            "di sulfide": "[#16X2H0][#16X2H0]",
            "epoxide": "[#6]1-[#6]-[#8]-1",
            "enamine": "[NX3][CX3]=[CX3]",
            "enol": "[OX2H][#6X3]=[#6]",
            "ester": "[#6][CX3](=O)[OX2H0][#6]",
            "ether": "[OD2]([#6])[#6]",
            "fluorine": "[F]",
            "furan": "[#6]1:[#6]:[#6]:[#6]:[#8]:1",
            "hydrogen": "[H]",
            "hydrazine": "[NX3][NX3]",
            "hydrazone": ["[NX3][NX2]=[*]", '[#6]-[#7](-[#6])-[#7]=[#6](-[#6])-[#6]'],
            "hydroxyl": "[OX2H]",
            "hydroxyl in alcohol": "[#6][OX2H]",
            "hydroxyl in carboxylic acid": "[OX2H][CX3]=[OX1]",
            "isonitrile": "[CX1-]#[NX2+]",
            "imidazoline": "[#6]1=[#7]-[#6]-[#6]-[#7]-1",
            "imide": "[CX3](=[OX1])[NX3H][CX3](=[OX1])",
            "imine": "[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]",
            "iminium": "[NX3+]=[CX3]",
            "iodine": "[#53]",
            "isoxazole": "[#6]1:[#6]:[#6]:[#7]:[#8]:1",
            "ketone": "[CX3]=[OX1]",
            "lactam": "[#8]=[#6]1-[#7](-[#6])-[#6]-[#6]-1",
            "peroxide": "[OX2,OX1-][OX2,OX1-]",
            "phenol": "[OX2H][cX3]:[c]",
            "phosphoric acid": "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
            "phosphoric ester": "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
            "primary alcohol": "[OX2H]",
            "primary amine": "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",
            "propane": "[#6]-[#6]-[#6]",
            "proton": "[H+]",
            "pyrazine": "[#6]1:[#6]:[#7]:[#6]:[#6]:[#7]:1",
            "quinazoline": "[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#6]:[#7]:2",
            "quinoline": "[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#6]:[#7]:2",
            "mono sulfide": "[#16X2H0][!#16]",
            "morpholine": "[#6]1-[#6]-[#7]-[#6]-[#6]-[#8]-1",
            "nitrate": "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]",
            "nitrile": "[NX1]#[CX2]",
            "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
            "nitroso": "[NX2]=[OX1]",
            "n-oxide": "[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
            "oxazolidinone": "[#8]=[#6]1-[#8]-[#6]-[#6]-[#7]-1",
            "secondary amine": "[NX3;H2,H1;!$(NC=O)]",
            "sulfate": "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
            "sulfamate": "[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]",
            "sulfamic acid": "[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]",
            "sulfenic acid": "[#16X2][OX2H,OX1H0-]",
            "sulfenate": "[#16X2][OX2H0]",
            "sulfide": "[#16X2H0]",
            "sulfonate": "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]",
            "sulfinate": "[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]",
            "sulfinic acid": "[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]",
            "sulfonamide": "[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]",
            "sulfone": "[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]",
            "sulfonic acid": "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
            "sulfoxide": "[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]",
            "sulfur": "[#16!H0]",
            "sulfuric acid ester": "[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]",
            "sulfuric acid diester": "[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]",
            "thioaldehyde": "[#16]=[#6H]-[#6]",
            "thioketone": "[#16]=[#6](-[#6])-[#6]", 
            "thioamide": "[NX3][CX3]=[SX1]",
            "thiol": "[#16X2H]",
            "vinylic carbon": "[$([CX3]=[CX3])]",
        }

        return functional_groups_smarts

    #------------------------- Property Declaration for GlobalChem ---------------------------#

    amino_acid_side_chains = property(_get_amino_acids)
    functional_groups_smiles = property(_get_functional_groups_smiles)
    functional_groups_smarts = property(_get_functional_groups_smarts)
    common_regex_patterns = property(_get_common_regex_patterns)
