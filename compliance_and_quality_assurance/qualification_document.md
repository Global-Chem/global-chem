Installation and Operational Qualification Testing for GlobalChem
=================================================================

This document will provide as a qualification document for `GlobalChem` being tested on multiple platforms and it's insurance
to it's dedicated functionality. 

## Table of Contents

<details>
 <summary>Click to open TOC</summary>

- [Introduction](#introduction)
- [System Architecture](#system-architecture)
- [Test Architecture](#test-architecture)
- [Installation Architecture](#installation-architecture)
- [GlobalChem Installation Qualifications - System Tests](#installation-qualification)
- [GlobalChem Core Operational Qualifications - System Tests](#operational-qualification)
- [Summary of Findings](#summary-findings)

</details>

# Introduction

This document will provide as a qualification document for `GlobalChem` being tested on multiple platforms and it's insurance
to it's dedicated functionality. `GlobalChem` is entirely written in one computing language, python, and has no extra dependencies
on other software. 

# System Architecture

`GlobalChem` is supported in Windows, Linux, and MacOS with backwards compatiblity. 

| GlobalChem Version | Latest Platform Support                  | Python Versions Supported  |
|--------------------|------------------------------------------|----------------------------|
| 1.6.1.5            | Linux x86_64, Windows 11, MacOS Monterey | 3.7, 3.8, 3.9              |

# Test Architecture

`GlobalChem` has a simplistic test architecture where all modules are tested accurately with `PyTest` using `Github Actions`
to facilitate the workflow management. Tests are performed on every "Pull Request", "Merge", or "Commit" into the `master`
branch. 

`FOSSA` is an open source software managemenment tool to check compatibility between open source licenses. 

| Test Software      | Continuous Integration & Deployment   | Coverage Software          | Legal Testing Software     |
|--------------------|---------------------------------------|----------------------------|----------------------------|
| PyTest             | Github Actions                        | Nose                       | FOSSA                      |


# Installation Architecture

`GlobalChem` is distributed on PyPi as a gzipped tar file and can be "pip", python installation service, installed directly
from a command line:

```python

pip install global-chem

```

| Distribution       | Version        | Version Hash (SHA256)                                            | File Type    |
|--------------------|----------------|------------------------------------------------------------------|--------------|
| PyPi               | 1.6.1.5        | 09eec09c0d6f2e8d8362662d8b559a966e8e409effe50f61e838a1a2401a1293 | Tar GZIP     |

# Installation Qualification

`GlobalChem` is tested and installed on `Windows`, `Linux` and `MacOS`.

```python

python setup.py install

```

`Linux` with Python `3.7`:

```bash

creating build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/PKG-INFO -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/SOURCES.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/dependency_links.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/not-zip-safe -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/requires.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/top_level.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
creating dist
creating 'dist/global_chem-1.6.1.5-py3.7.egg' and adding 'build/bdist.linux-x86_64/egg' to it
removing 'build/bdist.linux-x86_64/egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.7.egg
creating /opt/hostedtoolcache/Python/3.7.13/x64/lib/python3.7/site-packages/global_chem-1.6.1.5-py3.7.egg
Extracting global_chem-1.6.1.5-py3.7.egg to /opt/hostedtoolcache/Python/3.7.13/x64/lib/python3.7/site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed /opt/hostedtoolcache/Python/3.7.13/x64/lib/python3.7/site-packages/global_chem-1.6.1.5-py3.7.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

`Linux` with Python `3.8`:

```bash

creating build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/PKG-INFO -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/SOURCES.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/dependency_links.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/not-zip-safe -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/requires.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/top_level.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
creating dist
creating 'dist/global_chem-1.6.1.5-py3.8.egg' and adding 'build/bdist.linux-x86_64/egg' to it
removing 'build/bdist.linux-x86_64/egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.8.egg
creating /opt/hostedtoolcache/Python/3.8.12/x64/lib/python3.8/site-packages/global_chem-1.6.1.5-py3.8.egg
Extracting global_chem-1.6.1.5-py3.8.egg to /opt/hostedtoolcache/Python/3.8.12/x64/lib/python3.8/site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed /opt/hostedtoolcache/Python/3.8.12/x64/lib/python3.8/site-packages/global_chem-1.6.1.5-py3.8.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

`Linux` with Python `3.9`:

```bash

creating build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/PKG-INFO -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/SOURCES.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/dependency_links.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/not-zip-safe -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/requires.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
copying global_chem.egg-info/top_level.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
creating dist
creating 'dist/global_chem-1.6.1.5-py3.9.egg' and adding 'build/bdist.linux-x86_64/egg' to it
removing 'build/bdist.linux-x86_64/egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.9.egg
creating /opt/hostedtoolcache/Python/3.9.12/x64/lib/python3.9/site-packages/global_chem-1.6.1.5-py3.9.egg
Extracting global_chem-1.6.1.5-py3.9.egg to /opt/hostedtoolcache/Python/3.9.12/x64/lib/python3.9/site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed /opt/hostedtoolcache/Python/3.9.12/x64/lib/python3.9/site-packages/global_chem-1.6.1.5-py3.9.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

`Windows` with Python `3.7`

```bash

creating build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\PKG-INFO -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\SOURCES.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\dependency_links.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\not-zip-safe -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\requires.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\top_level.txt -> build\bdist.win-amd64\egg\EGG-INFO
creating dist
creating 'dist\global_chem-1.6.1.5-py3.7.egg' and adding 'build\bdist.win-amd64\egg' to it
removing 'build\bdist.win-amd64\egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.7.egg
creating c:\hostedtoolcache\windows\python\3.7.9\x64\lib\site-packages\global_chem-1.6.1.5-py3.7.egg
Extracting global_chem-1.6.1.5-py3.7.egg to c:\hostedtoolcache\windows\python\3.7.9\x64\lib\site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed c:\hostedtoolcache\windows\python\3.7.9\x64\lib\site-packages\global_chem-1.6.1.5-py3.7.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

`Windows` with Python `3.8`

```bash

creating build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\PKG-INFO -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\SOURCES.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\dependency_links.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\not-zip-safe -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\requires.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\top_level.txt -> build\bdist.win-amd64\egg\EGG-INFO
creating dist
creating 'dist\global_chem-1.6.1.5-py3.8.egg' and adding 'build\bdist.win-amd64\egg' to it
removing 'build\bdist.win-amd64\egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.8.egg
creating c:\hostedtoolcache\windows\python\3.8.10\x64\lib\site-packages\global_chem-1.6.1.5-py3.8.egg
Extracting global_chem-1.6.1.5-py3.8.egg to c:\hostedtoolcache\windows\python\3.8.10\x64\lib\site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed c:\hostedtoolcache\windows\python\3.8.10\x64\lib\site-packages\global_chem-1.6.1.5-py3.8.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

`Windows` with Python `3.9`

```bash

creating build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\PKG-INFO -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\SOURCES.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\dependency_links.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\not-zip-safe -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\requires.txt -> build\bdist.win-amd64\egg\EGG-INFO
copying global_chem.egg-info\top_level.txt -> build\bdist.win-amd64\egg\EGG-INFO
creating dist
creating 'dist\global_chem-1.6.1.5-py3.9.egg' and adding 'build\bdist.win-amd64\egg' to it
removing 'build\bdist.win-amd64\egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.9.egg
creating c:\hostedtoolcache\windows\python\3.9.12\x64\lib\site-packages\global_chem-1.6.1.5-py3.9.egg
Extracting global_chem-1.6.1.5-py3.9.egg to c:\hostedtoolcache\windows\python\3.9.12\x64\lib\site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed c:\hostedtoolcache\windows\python\3.9.12\x64\lib\site-packages\global_chem-1.6.1.5-py3.9.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5

```

`MacOS` with Python `3.7`


```bash

creating build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/PKG-INFO -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/SOURCES.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/dependency_links.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/not-zip-safe -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/requires.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/top_level.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
creating dist
creating 'dist/global_chem-1.6.1.5-py3.7.egg' and adding 'build/bdist.macosx-10.15-x86_64/egg' to it
removing 'build/bdist.macosx-10.15-x86_64/egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.7.egg
creating /Users/runner/hostedtoolcache/Python/3.7.13/x64/lib/python3.7/site-packages/global_chem-1.6.1.5-py3.7.egg
Extracting global_chem-1.6.1.5-py3.7.egg to /Users/runner/hostedtoolcache/Python/3.7.13/x64/lib/python3.7/site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed /Users/runner/hostedtoolcache/Python/3.7.13/x64/lib/python3.7/site-packages/global_chem-1.6.1.5-py3.7.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

`MacOS` with Python `3.8`

```bash

creating build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/PKG-INFO -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/SOURCES.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/dependency_links.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/not-zip-safe -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/requires.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/top_level.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
creating dist
creating 'dist/global_chem-1.6.1.5-py3.8.egg' and adding 'build/bdist.macosx-10.15-x86_64/egg' to it
removing 'build/bdist.macosx-10.15-x86_64/egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.8.egg
creating /Users/runner/hostedtoolcache/Python/3.8.12/x64/lib/python3.8/site-packages/global_chem-1.6.1.5-py3.8.egg
Extracting global_chem-1.6.1.5-py3.8.egg to /Users/runner/hostedtoolcache/Python/3.8.12/x64/lib/python3.8/site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed /Users/runner/hostedtoolcache/Python/3.8.12/x64/lib/python3.8/site-packages/global_chem-1.6.1.5-py3.8.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```
`MacOS` with Python `3.9`

```bash

creating build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/PKG-INFO -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/SOURCES.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/dependency_links.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/not-zip-safe -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/requires.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
copying global_chem.egg-info/top_level.txt -> build/bdist.macosx-10.15-x86_64/egg/EGG-INFO
creating dist
creating 'dist/global_chem-1.6.1.5-py3.9.egg' and adding 'build/bdist.macosx-10.15-x86_64/egg' to it
removing 'build/bdist.macosx-10.15-x86_64/egg' (and everything under it)
Processing global_chem-1.6.1.5-py3.9.egg
creating /Users/runner/hostedtoolcache/Python/3.9.12/x64/lib/python3.9/site-packages/global_chem-1.6.1.5-py3.9.egg
Extracting global_chem-1.6.1.5-py3.9.egg to /Users/runner/hostedtoolcache/Python/3.9.12/x64/lib/python3.9/site-packages
Adding global-chem 1.6.1.5 to easy-install.pth file

Installed /Users/runner/hostedtoolcache/Python/3.9.12/x64/lib/python3.9/site-packages/global_chem-1.6.1.5-py3.9.egg
Processing dependencies for global-chem==1.6.1.5
Finished processing dependencies for global-chem==1.6.1.5
```

# Operational Qualification

`GlobalChem` has a series of tests that will be performed by Github Actions bot on the PyTest Modules the developers have 
wrote. The tests are applied across `Linux`, `Windows`, and `MacOS` and can be performed with the following command:

```python

python -m nose --verbose --with-coverage -s -w tests/

```

Test Results Across all systems

```bash

Name                                                                                      Stmts   Miss  Cover
-------------------------------------------------------------------------------------------------------------
GlobalChem/__init__.py                                                                       33      0   100%
global_chem/environment/__init__.py                                                           0      0   100%
global_chem/environment/emerging_perfluoroalkyls.py                                          11      0   100%
global_chem/formulation/__init__.py                                                           0      0   100%
global_chem/formulation/excipients/__init__.py                                                0      0   100%
global_chem/formulation/excipients/biopharmaceutics_class_three/__init__.py                   0      0   100%
global_chem/formulation/excipients/biopharmaceutics_class_three/cimetidine_acyclovir.py      11      0   100%
global_chem/formulation/excipients/monoclonal_antibodies/__init__.py                          0      0   100%
global_chem/formulation/excipients/monoclonal_antibodies/monoclonal_antibodies.py            11      0   100%
global_chem/global_chem.py                                                                  417     80    81%
global_chem/interstellar_space/__init__.py                                                    0      0   100%
global_chem/interstellar_space/interstellar_space.py                                         11      0   100%
global_chem/materials/__init__.py                                                             0      0   100%
global_chem/materials/clay/__init__.py                                                        0      0   100%
global_chem/materials/clay/montmorillonite_adsorption.py                                     11      0   100%
global_chem/materials/polymers/__init__.py                                                    2      0   100%
global_chem/materials/polymers/common_monomer_repeating_units.py                             11      0   100%
global_chem/medicinal_chemistry/__init__.py                                                   0      0   100%
```

# Summary Findings

`GlobalChem` is qualified to install and perform on multiple operating systems with multiple python versions:

| GlobalChem Version | Platform \| Python Version | Installation Results | Operational Results |
|--------------------|----------------------------|----------------------|---------------------|
| 1.6.1.5            | Linux x86 \| 3.7           | PASS                 | PASS                |
| 1.6.1.5            | Linux x86 \| 3.8           | PASS                 | PASS                |
| 1.6.1.5            | Linux x86 \| 3.9           | PASS                 | PASS                |
| 1.6.1.5            | Windows 11 \| 3.7          | PASS                 | PASS                |
| 1.6.1.5            | Windows 11 \| 3.8          | PASS                 | PASS                |
| 1.6.1.5            | Windows 11 \| 3.9          | PASS                 | PASS                |
| 1.6.1.5            | Mac OS Monterey \| 3.7     | PASS                 | PASS                |
| 1.6.1.5            | Mac OS Monterey \| 3.8     | PASS                 | PASS                |
| 1.6.1.5            | Mac OS Monterey \| 3.9     | PASS                 | PASS                |

Currently, version `1.6.1.5` can perform adequately within the guidelines of 21 Part CFR 11 Compliance