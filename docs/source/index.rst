.. Global Chem documentation master file, created by
   sphinx-quickstart on Wed Aug 28 23:41:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GlobalChem's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. GlobalChem documentation master file, created by sphinx-quickstart on Tue Mar 24 16:12:38 2015.

GlobalChem
==========

.. sectionauthor:: Suliman Sharif <sharifsuliman1@gmail.com>

GlobalChem is a content variable store for all things chemistry.
Essentially it is a class object used to store attributes and return information for commonly used variables used across
multiple scripts.

>>> from global_chem import GlobalChem
>>> global_chem = GlobalChem()
>>> global_chem.amino_acid_side_chains
>>> global_chem.functional_groups_smiles

User guide
----------

A step-by-step guide to getting started with Cocktail Shaker.

.. toctree::
    :maxdepth: 2

    guide/installation
    guide/quickstart
    guide/contributing

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method.

.. toctree::
    :maxdepth: 2

    guide/content_store
    guide/functional_groups
    guide/amino_acids

Global Chem's license
---------------------

GlobalChem is released under the Mozilla Public License 2.0. This is a short, permissive software license that allows commercial use,
modifications, distribution, sublicensing and private use. Basically, you can do whatever you want with GlobalChem as long as
you include the original copyright and license in any copies or derivative projects.

.. _`MPL license`: https://github.com/Sulstice/global-chem/blob/master/LICENSE
