===============================================================
EXPOsan: EXPOsition of sanitation and resource recovery systems
===============================================================

.. image:: https://img.shields.io/pypi/l/exposan?color=blue&logo=UIUC&style=flat
   :target: https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
.. image:: https://img.shields.io/pypi/pyversions/exposan?style=flat
   :target: https://pypi.python.org/pypi/exposan
.. image:: https://img.shields.io/pypi/v/exposan?style=flat&color=blue
   :target: https://pypi.org/project/exposan/
.. image:: https://img.shields.io/appveyor/build/yalinli2/EXPOsan/main?label=build&logo=appveyor
   :target: https://github.com/QSD-Group/EXPOsan/tree/main
.. image:: https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg
   :target: https://qsdsan.readthedocs.io/en/latest/CODE_OF_CONDUCT.html


What is ``EXPOsan``?
--------------------
The ``EXPOsan`` package contains systems developed using the `QSDsan <https://github.com/QSD-Group/QSDsan>`_ package. Please refer to the individual README files within each system module for loading instructions and exemplary results.


Installation
------------
The easiest way is through ``pip``, in command/Anaconda prompt (Windows) or terminal (macOS):

.. code:: bash

    pip install exposan

If you need to update:

.. code::

    pip install -U exposan

or for a specific version (replace X.X.X with the version number):

.. code::

    pip install exposan==X.X.X

If you want to install the latest GitHub version at the `main branch <https://github.com/qsd-group/exposan>`_:

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

    pip install –r requirements.txt


.. note::
   
   Using the ``depth`` flag will only clone the main branch by default, if there are other branches and you want them, then use the ``no-single-branch`` flag, i.e.

   .. code:: bash

       git clone https://github.com/<YOUR_USERNAME>/EXPOsan.git --depth=1 --no-single-branch


Authors and Contributing
------------------------
``EXPOsan`` is developed by the Quantitative Sustainable Design Group led by members of the `Guest Group <http://engineeringforsustainability.com/>`_ and `other developers <https://github.com/QSD-Group/EXPOsan/graphs/contributors>`_ that have contributed to the repository. Yalin Li (`@yalinli2 <https://github.com/yalinli2>`_; `email <zoe.yalin.li@gmail.com>`_) is the currently maintainer of this repository and package.

Please refer to the `Contributing to QSDsan <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation for instructions and guidelines.


License information
-------------------
Please refer to the ``LICENSE.txt`` for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.