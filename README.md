About AtomCalc
==============
AtomCalc calculates the population evolution of a general n-level system in the presence of laser interaction. 
A system is defined by the levels, the lasers, and decay paths. For each of them exists one class that owns the corresponding properties.
The time evolution of the population of each level is then calculated with the 'fidelity' function that uses a Lindblad master equation approach.

This project is supposed to be expanded and should be seen as a construction fundament.
The methods used in the code are explained in my master thesis which can be obtained via the 5th institute of physics of the university of Stuttgart.


Installation
============
AtomCalc can be installed with ``pip install atomcalc``. Its dependencies are Matplotlib, SciPy, NumPy and QuTiP.