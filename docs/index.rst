.. AtomCalc documentation master file, created by
   sphinx-quickstart on Fri May  5 15:54:24 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AtomCalc's documentation!
====================================

AtomCalc calculates the population evolution of a general n-level system in the presence of laser interaction. 
A system is defined by the levels, the lasers, and decay paths. For each of them exists one class that owns the corresponding properties.
The time evolution of the population of each level is then calculated with the ``fidelity`` function that uses a Lindblad master equation approach.

This project is supposed to be expanded and should be seen as a construction fundament.
The methods used in the code are explained in my master thesis which can be obtained via the 5th institute of physics of the university of Stuttgart.

The code is documented in the `Code`_ section of this documentation. 
The whole python file can be accessed on the corresponding `GitHub page`_.
The tutorials can also be accessed on `GitHub`_ as Jupyter notebooks.



.. _GitHub page: https://github.com/AtomCalc/AtomCalc/tree/main/atomcalc
.. _GitHub: https://github.com/AtomCalc/AtomCalc/tree/main/docs/notebooks
.. _Code: https://atomcalc.github.io/AtomCalc/modules.html

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
   :caption: Code:

   modules

