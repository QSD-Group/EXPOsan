===============================================================
EXPOsan: EXPOsition of sanitation and resource recovery systems
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

.. GitHub test of the main branch
.. image:: https://github.com/QSD-Group/EXPOsan/actions/workflows/build-latest.yml/badge.svg?branch=main
   :target: https://github.com/QSD-Group/EXPOsan/actions/workflows/build-latest.yml

.. GitHub test of the stable branch (i.e., PyPI released package)
.. image:: https://github.com/QSD-Group/EXPOsan/actions/workflows/build-stable.yml/badge.svg?branch=stable
   :target: https://github.com/QSD-Group/EXPOsan/actions/workflows/build-stable.yml

.. Code of Conduct
.. image:: https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg
   :target: https://qsdsan.readthedocs.io/en/latest/CODE_OF_CONDUCT.html

.. AppVeyor test of the stable branch, not in active use
..
    image:: https://img.shields.io/appveyor/build/yalinli2/EXPOsan/main?label=build-stable&logo=appveyor
   :target: https://github.com/QSD-Group/EXPOsan/tree/stable


What is ``EXPOsan``?
--------------------
The ``EXPOsan`` package contains systems developed using the `QSDsan <https://github.com/QSD-Group/QSDsan>`_ package. Please refer to the individual README files within each system module for loading instructions and exemplary results.


Installation
------------
The easiest way is through ``pip``, in your command-line interface (e.g., Anaconda prompt, terminal):

.. code:: bash

    pip install exposan

If you need to upgrade:

.. code::

    pip install -U exposan

or for a specific version (replace X.X.X with the version number):

.. code::

    pip install exposan==X.X.X

If you want to install the latest GitHub version at the `main branch <https://github.com/qsd-group/exposan>`_ (note that you can still use the ``-U`` flag for upgrading):

.. code::

    pip install git+https://github.com/QSD-Group/EXPOsan.git


.. note::

   If this doesn't give you the newest ``exposan``, try ``pip uninstall exposan`` first.

   Also, you may need to update some ``exposan``'s dependency package (e.g., ' ``qsdsan``, ``biosteam``, ``thermosteam``) versions in order for the new ``exposan`` to run.


or other fork and/or branch (replace ``<USERNAME_OF_THE_FORK>`` and ``<BRANCH_NAME>`` with the desired fork and branch names)

.. code::

    pip install git+https://github.com/<USERNAME_OF_THE_FORK>/EXPOsan.git@<BRANCH_NAME>


You can also download the package from `PyPI <https://pypi.org/project/exposan/>`_.


To get the git version (use the ``depth`` flag to choose how many commit histories you want to clone):

.. code:: bash

    git clone https://github.com/QSD-Group/EXPOsan.git --depth=1

then navigate into the repository (``cd EXPOsan``) and install the required packages:

.. code:: bash

    pip install –r requirements.txt


.. note::
   
   Using the ``depth`` flag will only clone the main branch by default, if there are other branches and you want them, then use the ``no-single-branch`` flag, i.e.

   .. code:: bash

       git clone https://github.com/<YOUR_USERNAME>/EXPOsan.git --depth=1 --no-single-branch


Authors and Contributing
------------------------
``EXPOsan`` is developed by the Quantitative Sustainable Design Group led by members of the `Guest Group <http://engineeringforsustainability.com/>`_ and `other developers <https://github.com/QSD-Group/EXPOsan/graphs/contributors>`_ that have contributed to the repository. Yalin (`webpage <https://qsdsan.readthedocs.io/en/latest/authors/Yalin_Li.html>`_; `email <zoe.yalin.li@gmail.com>`_) is the currently maintainer of this repository and package.

Please refer to the `Contributing to QSDsan <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation for instructions and guidelines.


Stay Connected
--------------
If you would like to receive exciting news related to the QSDsan platform, you can subscribe to email updates using `this form <https://groups.webservices.illinois.edu/subscribe/154591>`_ (don't worry, you will be able to unsubscribe :)). Thank you in advance for your interest!


License information
-------------------
Please refer to the ``LICENSE.txt`` for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.