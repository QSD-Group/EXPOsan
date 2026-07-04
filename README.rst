===============================================================
EXPOsan: EXPOsition of Sanitation and Resource Recovery Systems
===============================================================

.. License
.. image:: https://img.shields.io/pypi/l/exposan?color=blue&logo=UIUC&style=flat
   :target: https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt

.. Tested Python version
.. image:: https://img.shields.io/pypi/pyversions/exposan?style=flat
   :target: https://pypi.python.org/pypi/exposan

.. PyPI version
.. image:: https://img.shields.io/pypi/v/exposan?style=flat&color=blue
   :target: https://pypi.org/project/exposan

.. Zenodo release DOI
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.20256578.svg
   :target: https://doi.org/10.5281/zenodo.20256578

.. GitHub test of the main branch
.. image:: https://github.com/QSD-Group/EXPOsan/actions/workflows/build-latest.yml/badge.svg?branch=main
   :target: https://github.com/QSD-Group/EXPOsan/actions/workflows/build-latest.yml

.. QSDsan documentation
.. image:: https://readthedocs.org/projects/qsdsan/badge/?version=latest
   :target: https://qsdsan.readthedocs.io/en/latest


What is ``EXPOsan``?
--------------------
The ``EXPOsan`` package contains systems developed using the `QSDsan <https://github.com/QSD-Group/QSDsan>`_ package. A list of the developed systems is included in QSDsan's `documentation <https://qsdsan.readthedocs.io/en/latest/Developed_Systems.html>`_. Please refer to the individual README files within each system module for loading instructions and exemplary results.


Installation
------------
``EXPOsan`` requires Python 3.12 or newer. The easiest way to install ``EXPOsan`` is through ``pip`` in a command-line interface (e.g., terminal, PowerShell, etc.):

.. code:: bash

    pip install exposan

This installs the core package. A few optional features (e.g., routing optimization
in the ``hap`` module) need extra packages.
You can install complete dependencies with:

.. code::

    pip install exposan[complete]

or install the specific package named in the error message when you first hit
a missing feature.

If you need to upgrade:

.. code::

    pip install -U exposan

or for a specific version (replace X.X.X with the version number):

.. code::

    pip install exposan==X.X.X

If you want to install the latest GitHub version at the `main branch <https://github.com/QSD-Group/EXPOsan>`_ (note that you can still use the ``-U`` flag for upgrading):

.. code::

    pip install git+https://github.com/QSD-Group/EXPOsan.git


or other fork and/or branch (replace ``<USERNAME_OF_THE_FORK>`` and ``<BRANCH_NAME>`` with the desired fork and branch names)

.. code::

    pip install git+https://github.com/<USERNAME_OF_THE_FORK>/EXPOsan.git@<BRANCH_NAME>


You can also download the package from `PyPI <https://pypi.org/project/exposan/>`_.


To get the git version (use the ``depth`` flag to choose how many commit histories you want to clone):

.. code:: bash

    git clone https://github.com/QSD-Group/EXPOsan.git --depth=1

then navigate into the repository (``cd EXPOsan``) and install the required packages:

.. code:: bash

    pip install -e ".[dev]"


.. note::
   
   Using the ``depth`` flag will only clone the main branch by default, if there are other branches and you want them, then use the ``no-single-branch`` flag, i.e.

   .. code:: bash

       git clone https://github.com/<YOUR_USERNAME>/EXPOsan.git --depth=1 --no-single-branch


Authors and Contributing
------------------------
``EXPOsan`` and ``QSDsan`` are developed by the Quantitative Sustainable Design Group and the broader community. `Yalin Li <https://github.com/yalinli2>`_ is the currently maintainer. See `commit history <https://github.com/QSD-Group/EXPOsan/graphs/contributors>`_ for contributors who have contributed to the repository. 

If you want to contribute to ``EXPOsan``, please refer to the `Contributing Guidelines <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation for instructions and guidelines.


License information
-------------------
Please refer to the ``LICENSE.txt`` for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.
