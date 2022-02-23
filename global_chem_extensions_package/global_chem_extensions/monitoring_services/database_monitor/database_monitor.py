#!/usr/bin/env python3
#
# GlobalChemExtensions - Database Monitor
#
# ----------------------------------------

# Imports
# -------

import urllib.request

class DatabaseMonitor(object):

    __version__ = '0.0.1'

    def __init__(self):

        self.urls, self.names = self.prepare_urls()

    def prepare_urls(self):

        '''

        Prepare the URL

        Seems meaningless but can expand this function later on.

        '''

        contents = urllib.request.urlopen("https://raw.githubusercontent.com/Sulstice/Uptime-Cheminformatics/master/.upptimerc.yml").readlines()
        contents = [ i.decode('utf-8') for i in contents ]

        urls = []
        names = []

        for i in contents:
            if 'name' in i and \
                '#' not in i and 'chemistrydb.com' not in i \
                and 'Cheminformatic Database Statuses' not in i:
                i = i.strip('\n').split('name:')[1]
                names.append(i)

            if 'url' in i:
                i = i.strip('\n').split('url:')[1]
                urls.append(i)

        return urls, names

    def heartbeat(self):

        '''

        Perform the heartbeat check

        '''

        successes = {}
        failures = {}

        for i in range(0, len(self.urls)):

            url = self.urls[i]
            name = self.names[i]
            print (name)

            try:
                response = urllib.request.urlopen(url)
                status_code = response.getcode()

                if '200' == status_code:
                    successes[name] = 'Up'
            except:
                failures[name] = 'Down'

        return successes, failures


