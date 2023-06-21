# %%
# Imports
import numpy as np
import qutip
import matplotlib.pylab as plt
import scipy


def plot_pulse(pulse, tlist):
    """
    A function to plot a pulse function.

    Args:
        pulse (function): value for :attr:`pulse`
        tlist (list): value for :attr:`tlist`

    Attributes:
        pulse (function): A time dependent pulse function.
        tlist (list): A list of points in time where the pulse should be plotted.

    """
    fig, ax = plt.subplots()
    pulse = np.array([pulse(t) for t in tlist])
    tlist = np.array(tlist)
    ax.set_yticks([0, 1], labels=[r"0", r"$\Omega$"])
    ax.set_ylim([-0.1, 1.1])
    ax.plot(tlist, pulse / max(pulse), label=r"$\Omega(t)$")
    ax.legend()
    ax.set_xlabel(r"Time")
    ax.set_ylabel(r"Pulse amplitude")
    plt.grid(linestyle=(0, (5, 10)), axis="both")
    plt.show(fig)
