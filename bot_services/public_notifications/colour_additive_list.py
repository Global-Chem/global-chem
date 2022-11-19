#!/usr/bin/env python3
#
# GlobalChem - Bot Colour Additive List
#
# -------------------------------------

# Imports
# -------

import os
import json
import requests
import pandas as pd

from bs4 import BeautifulSoup

class BotColourAdditiveList(object):

    def __init__(self):

        self.name = 'bot_colour_additive_list'
        self.master_url = "https://www.fda.gov/industry/color-additive-inventories/color-additive-status-list"

        self.fda_lists = []
        self.global_chem_lists = []

    def get_fda_reported_lists(self):

        '''
        Get the FDA Lists of Each Seven Colour Additives.
        '''

        page = requests.get(self.master_url)
        soup = BeautifulSoup(page.content, "html.parser")

        fda_lists = soup.find_all("ul", class_="ref")

        fda_current_status = []

        for fda_list in fda_lists:

            molecules = []
            items = fda_list.text.split('\n')
            items = [i for i in items if i]

            for item in items:

                molecule = item.split('-')[0].lower().strip()

                if 'c.i.' in molecule and not 'orange' in molecule:
                    molecule = molecule.split(' ')[1]

                import re
                molecule = re.sub(r'[^\x00-\x7f]',r'', molecule)

                if molecule == 'beta' or molecule == 'beta-carotene':
                    molecule = "beta-apo-8'-carotenal"

                if molecule == 'chromium':
                    molecule = 'chromium-cobalt-aluminum oxide'

                if molecule == 'fd&c lakes (see preface of this exhibit)':
                    molecule = 'fd&c lakes'

                if molecule == 'd&c lakes (see preface of this exhibit)':
                    molecule = 'd&c lakes'

                if molecule == 'd&c brown #1 only external cosmetics':
                    molecule = 'd&c brown #1'

                if 'ext.' in molecule:
                    continue

                if molecule == 'beet powder (dehydrated beets)':
                    molecule = 'beet powder'

                if molecule == 'beet juice (as vegetable juice)':
                    molecule = 'beet juice'

                if molecule == 'tagetes (aztec marigold) meal and extract':
                    molecule = 'tagetes'

                if molecule == 'beta carotene, natural and synthetic':
                    molecule = 'beta carotene'

                if molecule == 'grape skin extract (enocianina)':
                    molecule = 'grape skin extract'

                if molecule == 'paprika & paprika oleoresin':
                    molecule = 'paprika'

                if molecule == 'lycopene, tomato extract or concentrate':
                    molecule = 'lycopene'

                if molecule == 'turmeric & turmeric oleoresin':
                    molecule = 'turmeric'

                if molecule == 'cottonseed flour, toasted partially defatted cooked':
                    molecule = 'cottonseed flour'

                if molecule == 'mica-based pearlescent pigment':
                    molecule = 'mica-based pearlescent pigment'

                if molecule == 'mica':
                    molecule = 'mica'

                if 'spirulina extract' in molecule:
                    molecule = 'spirulina extract'

                if molecule == 'guanine (pearl essence)':
                    molecule = 'guanine'

                if molecule == 'ferric ammonium ferrocyanide (iron blue)':
                    molecule = 'ferric ammonium ferrocyanide'

                if molecule == 'potassium sodium copper chlorophyllin (chlorophyllin copper complex)':
                    molecule = 'potassium sodium copper chlorophyllin'

                if molecule == 'ferric ferrocyanide (iron blue)':
                    molecule = 'ferric ferrocyanide'

                if molecule == 'ultramarines (blue, green, pink, red & violet)':
                    molecule = 'ultramarines'

                if molecule == 'guaiazulene (azulene)':
                    molecule = 'guaiazulene'

                if molecule == 'disodium edta':
                    molecule = 'disodium edta-copper'

                if 'd&c yellow # 8' in molecule:
                    molecule = 'd&c yellow # 8'

                if molecule == 'd&c black # 4':
                    molecule = 'd&c black #4'

                if molecule == 'vinyl alcohol/methyl methacrylate':
                    molecule = 'vinyl alcohol'

                if molecule == 'chlorophyllin':
                    molecule = 'chlorophyllin-copper complex'

                molecules.append(molecule)

            fda_current_status.append(molecules)

        if len(fda_current_status) != 7:

            print ("Global-Chem Notification: Fda List Missing")

            raise ValueError

        self.fda_lists = fda_current_status

    def get_global_chem_lists(self):

        '''

        Check the Global-Chem Resource

        '''

        data = pd.read_csv(
            'https://raw.githubusercontent.com/Sulstice/global-chem/development/global_chem/global_chem_outputs/global_chem.tsv',
            sep='\t',
            names=['Name', 'SMILES', 'Node', 'Category', 'Node Path'],
            skiprows=1
        )

        list_one = []
        list_two = []
        list_three = []
        list_four = []
        list_five = []
        list_six = []
        list_seven = []

        for row in data.iterrows():

            node_path = row[1]['Node Path']
            name = row[1]['Name'].strip().lower()

            if 'color_additives' in node_path:

                if 'fda_list_one' in node_path:
                    list_one.append(name)
                if 'fda_list_two' in node_path:
                    list_two.append(name)
                if 'fda_list_three' in node_path:
                    list_three.append(name)
                if 'fda_list_four' in node_path:
                    list_four.append(name)
                if 'fda_list_five' in node_path:
                    list_five.append(name)
                if 'fda_list_six' in node_path:
                    list_six.append(name)
                if 'fda_list_seven' in node_path:
                    list_seven.append(name)

        global_chem_lists = [
            list_one,
            list_two,
            list_three,
            list_four,
            list_five,
            list_six,
            list_seven
        ]

        self.global_chem_lists = global_chem_lists

    def check_list_status(self):

        '''
        File a Github Issue on the global chem Repo if there is a misalignment.
        '''

        differences = []

        for i, chemical_list in enumerate(self.fda_lists):

            # Skip the Last List for Now since that will require more data curation

            if i == 6:
                continue

            difference = set(chemical_list).symmetric_difference(set(self.global_chem_lists[i]))

            for j in difference:

                if 'beta' in j or 'mica' in j:
                    pass
                else:
                    differences.append(j)


        if len(differences) > 1:

            token = os.getenv('GITHUB_TOKEN')
            
            headers = {
                "Accept": 'Accept: application/vnd.github+json',
                "Authorization" : "token {}".format(token)
            }

            data = {
                "title": "Public Notification System: Colour Additive List Changes",
                "body": "Differences found that need to be inspected: %s" % differences

            }

            username = 'Sulstice'
            repository_name = 'global-chem'

            url = "https://api.github.com/repos/{}/{}/issues".format(
                username,
                repository_name
            )

            response = requests.post(
                url,
                data=json.dumps(data),
                headers=headers
            )
            print (response.content)


if __name__ == '__main__':

    bot = BotColourAdditiveList()
    bot.get_fda_reported_lists()
    bot.get_global_chem_lists()
    bot.check_list_status()
