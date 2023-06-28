import numpy as np
import qutip
import matplotlib.pylab as plt
import scipy


class Level:
    """
    An object that describes an energy level.

    Args:
        energy (number): value for :attr:`energy`


    Attributes:
        energy (number): The energy of the level given as angular frequency.

    Example:

        >>> Level(5)
        This is a Level object with an energy of 5, given as an angular frequency.
    """

    def __init__(self, energy):
        self.energy = energy
