About AtomCalc
==============
[![Docs](https://github.com/AtomCalc/AtomCalc/actions/workflows/documentation.yaml/badge.svg)](https://atomcalc.github.io/AtomCalc/)

AtomCalc simulates the interaction between multi-level atoms and laser fields. It calculates how the population of the electronic levels of an atom changes if laser pulses are applied. To allow for realistic simulations, the motional state of the atom and decay channels can be taken into account.
A system is defined by the levels, the lasers, and decay paths. For each of them exists one class that owns the corresponding properties.
The time evolution of the population of each level is then calculated with the ``simulate`` function that uses a Lindblad master equation approach.

This project is supposed to be expanded and should be seen as a construction fundament.
The methods used in the code are explained in my master thesis which can be obtained via the 5th institute of physics of the university of Stuttgart.


Installation
============
AtomCalc can be installed with ``pip install atomcalc``. Its dependencies are Matplotlib, SciPy, NumPy and QuTiP.