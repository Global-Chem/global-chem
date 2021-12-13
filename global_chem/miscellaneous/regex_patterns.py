#!/usr/bin/env python3
#
# GlobalChem - Common Regex Patterns
#
# -----------------------------------

class CommonRegexPatterns(object):

    def __init__(self):

        pass

    def get_common_regex_patterns(self):

        regex_patterns = {
            'mol2': '^@<\w+?>\w+?\n[COMPOUND_ID]\n(.|\n)*?@<TRIPOS>SUBSTRUCTURE\n.*?\n'
        }

        return regex_patterns
