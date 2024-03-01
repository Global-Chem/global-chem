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

    def __init__(self, choice='chemistry'):

        self.choice = choice
        self.base_url = ''
        self._determine_parameters()
        self.urls, self.names = self.prepare_urls()

    def _determine_parameters(self):

        '''

        Determine the parameters that comes from the choice of chemistry or bacteria.

        '''

        if self.choice == 'chemistry':
            self.base_url = 'https://raw.githubusercontent.com/Sulstice/Uptime-Cheminformatics/master/.upptimerc.yml'
        elif self.choice == 'bacteria':
            self.base_url = 'https://raw.githubusercontent.com/Sulstice/Uptime-Bacteria/master/.upptimerc.yml'

    def prepare_urls(self):

        '''

        Prepare the URL

        Seems meaningless but can expand this function later on.

        '''

        contents = urllib.request.urlopen(self.base_url).readlines()
        contents = [ i.decode('utf-8') for i in contents ]

        urls = []
        names = []

        for i in contents:
            if 'name' in i and \
                '#' not in i and '.com' not in i \
                and ' Database Statuses' not in i:
                i = i.strip('\n').split('name:')[1]
                names.append(i)

            if 'url' in i:
                i = i.strip('\n').split('url:')[1]
                urls.append(i)

        return urls, names

    def heartbeat(self, verbose=False):

        '''

        Perform the heartbeat check

        Arguments:
            verbose (Boolean): if the user wants the verbose flag.

        '''

        successes = {}
        failures = {}

        for i in range(0, len(self.urls)):

            url = self.urls[i]
            name = self.names[i]
            status_code = ''

            try:
                response = urllib.request.urlopen(url)
                status_code = response.getcode()
                successes[name] = 'Up'
            except:
                failures[name] = 'Down'

            if verbose:
                print ("Server Name: %s, Status Code: %s" % (name, status_code))

        return successes, failures


