.. _contributing:

Contributing
============

.. sectionauthor:: Suliman Sharif <sharifsuliman1@gmail.com>

Contributing
------------

Contributions of any kind to global-chem are greatly appreciated! All contributions are welcome, no matter how big or small.

If you are able to contribute changes yourself, just fork the `source code`_ on GitHub, make changes, and file a pull
request.

Quick Guide to Contributing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. `Fork the GlobalChem repository on GitHub`_, then clone the fork to your local machine::

    git clone https://github.com/<username>/global-chem.git

2. Install the development requirements::

    cd global-chem
    pip install -r requirements/development.txt

3. Create a new branch for your changes::

    git checkout -b <name-for-changes>

4. Make your changes or additions. Ideally create some tests and ensure they pass.

5. Commit your changes and push to your fork on GitHub::

    git add .
    git commit -m "<description-of-changes>"
    git push origin <name-for-changes>

4. `Submit a pull request`_.

Tips
~~~~

- Follow the `PEP8`_ style guide.
- Include docstrings as described in `PEP257`_.
- Try and include tests that cover your changes.
- Try to write `good commit messages`_.
- Consider `squashing your commits`_ with rebase.
- Read the GitHub help page on `Using pull requests`_.

.. _`source code`: https://github.com/Sulstice/Cocktail-Shaker.git
.. _`Fork the CIRpy repository on GitHub`: https://github.com/Sulstice/Cocktail-Shaker/fork
.. _`Submit a pull request`: https://github.com/Sulstice/Cocktail-Shaker/compare/
.. _`squashing your commits`: http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html
.. _`PEP8`: https://www.python.org/dev/peps/pep-0008
.. _`PEP257`: https://www.python.org/dev/peps/pep-0257
.. _`good commit messages`: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
.. _`Using pull requests`: https://help.github.com/articles/using-pull-requests
