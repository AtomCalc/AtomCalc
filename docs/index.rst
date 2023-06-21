.. AtomCalc documentation master file, created by
   sphinx-quickstart on Fri May  5 15:54:24 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AtomCalc's documentation!
====================================

AtomCalc simulates the interaction between multi-level atoms and laser fields. It calculates how the population of the electronic levels of an atom changes if laser pulses are applied. To allow for realistic simulations, the motional state of the atom and decay channels can be taken into account.
A system is defined by the levels, the lasers, and decay paths. For each of them exists one class that owns the corresponding properties.
The time evolution of the population of each level is then calculated with the ``simulate`` function that uses a Lindblad master equation approach.

This project is supposed to be expanded and should be seen as a construction fundament.
The methods used in the code are explained in my master thesis which can be obtained via the 5th institute of physics of the university of Stuttgart.

The code is documented in the `Modules`_ section of this documentation. 
The whole python file can be accessed on the corresponding `GitHub page`_.
The tutorials can also be accessed on GitHub as `Jupyter notebooks`_.



.. _GitHub page: https://github.com/AtomCalc/AtomCalc
.. _Jupyter notebooks: https://github.com/AtomCalc/AtomCalc/tree/main/docs/notebooks
.. _Modules: https://atomcalc.github.io/AtomCalc/modules.html

Installation
============

AtomCalc can be installed with ``pip install atomcalc``. Its dependencies are Matplotlib, SciPy, NumPy and QuTiP.

.. toctree::
   :numbered:
   :maxdepth: 2
   :caption: Tutorials:

   notebooks/Example1
   notebooks/Example2
   notebooks/Example3

.. toctree::
   :maxdepth: 2
   :caption: Classes and Functions:

   modules

