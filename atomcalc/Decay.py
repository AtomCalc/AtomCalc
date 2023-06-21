import numpy as np
import qutip
import matplotlib.pylab as plt
import scipy


class Decay:
    """An object that describes the decay of the system with the decay rates and the respective transitions.

    Args:
        rates (list): value for :attr:`rates`
        final_states (list): value for :attr:`final_states`


    Attributes:
        rates (list): A list of decay rates.
        final_states (list): A list of tupels of :class:`Level` objects that assign the decay rates to a corresponding transition.

    Example:

        >>> Decay([0, 1], [[Level(20), Level(0)], [Level(5), Level(0)]])
        The transition between Level(20) and Level(0) is assigned a decay rate of 0. The transition between Level(5) and Level(0) is assigned a decay rate of 1.
    """

    def __init__(
        self, rates, final_states
    ):  # final_states is a Level-list with two-level-couples in each entry. rates is the respective decay-rate-list to each couple
        if type(rates) != list or type(final_states) != list:
            raise TypeError("rates and final_states need to be a list")
        self.rates = np.array(rates)
        self.final_states = np.array(final_states)
