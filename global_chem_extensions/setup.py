#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# package setup
#
# ------------------------------------------------

# imports
# -------
import os

# requirements
# ------------

with open('requirements.txt') as f:
    REQUIREMENTS = f.read().strip().split('\n')

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
    long_description = open('README.md').read()
else:
    long_description = 'GlobalChemExtensions - Your Extension functionality for GlobalChem'


# exec
# ----
setup(
    name="global_chem_extensions",
    version="1.0.2",
    packages=find_packages(),
    license='MPL 2.0',
    author="Suliman Sharif",
    author_email="sharifsuliman1@gmail.com",
    url="https://www.github.com/Sulstice/global-chem-extensions",
    install_requires=REQUIREMENTS,
    extras_require={
        'cheminformatics': [
            'partialsmiles', 'pysmiles', 'deepsmiles',
            'selfies', 'molvs', 'flask', 'plotly', 'kaleido',
            'bokeh', 'molpdf', 'dimorphite_dl',
            'scaffoldgraph'
        ],
        'bioinformatics': [
            'biopython', 'dna_features_viewer', 'biopandas',
            'pypdb'
        ],
        'development_operations': [''],
        'quantum_chemistry': ['moly', 'kaleido', 'pyyaml==3.13'],
        'forcefields': [
            'rdkit-pypi', 'partialsmiles', 'pysmiles', 'deepsmiles',
            'selfies', 'molvs'
        ],
        'graphing': ['plotly'],
    },
    long_description=long_description,
    long_description_content_type='text/markdown',
    zip_safe=False,
    keywords='smiles molecules chemistry rdkit plotly cheminformatics chemoinformatics organic python',
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
)
