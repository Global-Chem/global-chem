#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# package setup
#
# ------------------------------------------------

# imports
# -------
import os

# config
# ------
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

TEST_REQUIREMENTS = [
    'pytest',
    'pytest-runner'
]

if os.path.exists('README.md'):
    long_description = open('README.md', 'r', encoding='utf-8').read()
else:
    long_description = 'GlobalChem - Your Chemical Knowledge Graph for Chemistry'

# exec
# ----
setup(
    name="global_chem",
    version="1.6.1.3",
    packages=find_packages(),
    license='MPL 2.0',
    author="Suliman Sharif",
    author_email="sharifsuliman1@gmail.com",
    url="https://www.github.com/Sulstice/global-chem",
    install_requires=[],
    long_description=long_description,
    long_description_content_type='text/markdown',
    extras_require={
        'graphing': ['global-chem-extensions[graphing]'],
        'forcefields': ['global-chem-extensions[forcefields]'],
        'bioinformatics': ['global-chem-extensions[bioinformatics]'],
        'cheminformatics': ['global-chem-extensions[cheminformatics]'],
        'quantum_chemistry': ['global-chem-extensions[quantum_chemistry]'],
        'development_operations': ['global-chem-extensions[development_operations]'],
        'all': [
            'global-chem-extensions',
            'global-chem-extensions[graphing]',
            'global-chem-extensions[forcefields]',
            'global-chem-extensions[bioinformatics]',
            'global-chem-extensions[cheminformatics]',
            'global-chem-extensions[quantum_chemistry]',
            'global-chem-extensions[development_operations]',
        ]
    },
    zip_safe=False,
    keywords='smiles molecules chemistry organic iupac',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    tests_require=TEST_REQUIREMENTS,
)
