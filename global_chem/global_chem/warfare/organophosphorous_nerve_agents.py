#!/usr/bin/env python3
#
# GlobalChem - OrganoPhosphorous Nerve Agents
#
# -------------------------------------------

class OrganoPhosphorousNerveAgents(object):

    def __init__(self):

        self.name = 'organophosphorous_nerve_agents'

    @staticmethod
    def get_smiles():

        smiles = {
            'tabun': 'CCOP(=O)(C#N)N(C)C',
            'sarin': 'CC(C)OP(=O)(C)F',
            'soman': 'CC(C(C)(C)C)OP(=O)(C)F',
            'cyclosarin': 'CP(=O)(OC1CCCCC1)F',
            'vx': 'CCOP(C)(=O)SCCN(C(C)C)C(C)C',
            'russian vx': 'CCN(CC)CCSP(=O)(C)OCC(C)C',
            'mirzayanov-a230': 'CCN(CC)C(C)=N[P](C)(F)=O',
            'mirzayanov-a232': 'CCN(CC)C(\C)=N\P(F)(=O)OC',
            'mirzayanov-a234': r'CCOP(F)(=O)\N=C(/C)\N(CC)CC',
            'hoenig-a230': r'Cl/C(F)=N/OP(F)(OCCCl)=O',
            'hoenig-a232': r'Cl/C(F)=N/OP(F)(OC(C)CCl)=O',
            'hoenig-a234': r'Cl/C(F)=N/OP(F)(OC(C)C(C)Cl)=O',
            'novichok-5': 'FP1OC(C)CO1',
            'novichok-7': 'FP1OC(C)C(C)O1',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'tabun': '[#6]-[#6]-[#8]-[#15](=[#8])(-[#6]#[#7])-[#7](-[#6])-[#6]',
            'sarin': '[#6]-[#6](-[#6])-[#8]-[#15](=[#8])(-[#6])-[#9]',
            'soman': '[#6]-[#6](-[#6](-[#6])(-[#6])-[#6])-[#8]-[#15](=[#8])(-[#6])-[#9]',
            'cyclosarin': '[#6]-[#15](=[#8])(-[#8]-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1)-[#9]',
            'vx': '[#6]-[#6]-[#8]-[#15](-[#6])(=[#8])-[#16]-[#6]-[#6]-[#7](-[#6](-[#6])-[#6])-[#6](-[#6])-[#6]',
            'russian vx': '[#6]-[#6]-[#7](-[#6]-[#6])-[#6]-[#6]-[#16]-[#15](=[#8])(-[#6])-[#8]-[#6]-[#6](-[#6])-[#6]',
            'mirzayanov-a230': '[#6]-[#6]-[#7](-[#6]-[#6])-[#6](-[#6])=[#7]-[#15](-[#6])(-[#9])=[#8]',
            'mirzayanov-a232': '[#6]-[#6]-[#7](-[#6]-[#6])/[#6](-[#6])=[#7]/[#15](-[#9])(=[#8])-[#8]-[#6]',
            'mirzayanov-a234': '[#6]-[#6]-[#8]-[#15](-[#9])(=[#8])/[#7]=[#6](\[#6])-[#7](-[#6]-[#6])-[#6]-[#6]',
            'hoenig-a230': '[#17]/[#6](-[#9])=[#7]/[#8]-[#15](-[#9])(-[#8]-[#6]-[#6]-[#17])=[#8]',
            'hoenig-a232': '[#17]/[#6](-[#9])=[#7]/[#8]-[#15](-[#9])(-[#8]-[#6](-[#6])-[#6]-[#17])=[#8]',
            'hoenig-a234': '[#17]/[#6](-[#9])=[#7]/[#8]-[#15](-[#9])(-[#8]-[#6](-[#6])-[#6](-[#6])-[#17])=[#8]',
            'novichok-5': '[#9]-[#15]1-[#8]-[#6](-[#6])-[#6]-[#8]-1',
            'novichok-7': '[#9]-[#15]1-[#8]-[#6](-[#6])-[#6](-[#6])-[#8]-1',
        }

        return smarts
