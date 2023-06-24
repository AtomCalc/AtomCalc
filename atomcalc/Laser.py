import numpy as np
import qutip
import matplotlib.pylab as plt
import scipy


class Laser:
    """
    An object that describes values for a laser: Rabi frequency, detuning, and coupled states.
    One laser only affects the one specific transition that is defined in the couple parameter.

    Args:
        rabifreq (number): value for :attr:`rabifreq`
        detuning (number): value for :attr:`detuning`
        couple (list): value for :attr:`couple`
        pulse (None or function or list): value for :attr:`pulse`


    Attributes:
        rabifreq (number): Rabi frequency of the laser as angular frequency.
        detuning (number): The detuning of the laser as angular frequency.
        frequency (number): The to the frequency of the laser corresponding energy. Calculated by energy level difference + detuning.
        couple (list): A tupel of :class:`Level` objects that assigns the laser to this transition.
        pulse (None or function or list): A time dependent function of the Rabi frequency of the laser OR a list of numbers describing a Rabi frequency pulse.

    Example:

        >>> Laser(1, 100, [Level(0),Level(20)])
        The transition between Level(20) and Level(0) is assigned a laser with Rabi frequency of 1 and a frequency of 100.

    Note:
        Level couples need to be sorted by energy in ascending order.
    """

    def __init__(self, rabifreq, detuning, couple, polarization=None, pulse=None):
        self.rabi = rabifreq
        self.couple = couple
        self.detuning = detuning
        self.frequency = np.abs(couple[0].energy - couple[1].energy) + detuning
        self.polarization = polarization  # Not yet included. A list of a normalized E-field vector in the laser coordinate system, a theta_k and a theta_p (in degrees) and a pair [m_i,m_f].
        self.pulse = pulse
