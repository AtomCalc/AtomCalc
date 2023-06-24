import numpy as np
import qutip
import matplotlib.pylab as plt
import scipy


class System:
    """An object that inherits all parameters used for simulation of the system.

    Args:
        levels (list): value for :attr:`levels`
        lasers (list): value for :attr:`lasers`
        decay (class object): value for :attr:`decay`


    Attributes:
        levels (list): A list of :class:`Level` objects.
        lasers (list): A list of :class:`Laser` objects.
        decay (class object): A :class:`Decay` object.
        dim (number): Number of levels.

    Example:
        >>> level1 = Level(0)
        >>> level2 = Level(20)
        >>> level3 = Level(100)
        >>> laser1 = Laser(1, 120, [level1,level3])
        >>> laser2 = Laser(1, 100, [level2,level3])
        >>> decay = Decay([0],[[level3,level1]])
        >>> system = System([level1, level2, level3], [laser1,laser2], decay)

    Note:
        Levels need to be sorted by energy in ascending order.
    """

    def __init__(self, levels, lasers, decay):
        self.lasers = lasers
        self.levels = levels
        self.dim = len(levels)
        self.decay = decay

    def draw(self):
        # Levels
        for i in range(len(self.levels)):
            x = [0.01 + 0.2 * i, 0.19 + 0.2 * i]
            y = [self.levels[i].energy, self.levels[i].energy]
            plt.plot(x, y, linestyle="-", linewidth=1, label="Level {}".format(i + 1))
            plt.text(0.5 * (x[0] + x[1]), y[0] + 0.002, "{}".format(i + 1))
        plt.xlim([-0.25, 0.2 * len(self.levels) + 0.25])
        # Detuning
        detuning = [[]] * len(self.lasers)
        for i in range(len(self.lasers)):
            max_couple = max(
                self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy
            )
            min_couple = min(
                self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy
            )
            detuning[i] = min_couple + self.lasers[i].frequency - max_couple
        for i in range(len(self.lasers)):
            max_couple = max(
                self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy
            )
            if max_couple == self.lasers[i].couple[0].energy:
                i1 = self.levels.index(self.lasers[i].couple[0])
            else:
                i1 = self.levels.index(self.lasers[i].couple[1])
            x = [0.01 + 0.2 * i1, 0.19 + 0.2 * i1]
            y = [max_couple + detuning[i], max_couple + detuning[i]]
            plt.plot(
                x,
                y,
                linestyle="--",
                label="D {} = {:.2e}".format(i + 1, detuning[i]),
            )
        # Lasers
        for i in range(len(self.lasers)):
            max_couple = max(
                self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy
            )
            min_couple = min(
                self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy
            )
            i1 = self.levels.index(self.lasers[i].couple[0])
            i2 = self.levels.index(self.lasers[i].couple[1])
            x = [0.1 + 0.2 * i1, 0.1 + 0.2 * i2]
            y = [min_couple, max_couple + detuning[i]]
            plt.plot(x, y, linestyle="-", label="Laser {}".format(i + 1))
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.show()

    def _plot_population(self, result, dim):
        fig, ax = plt.subplots()
        for i in range(dim):
            if i != 12:
                ax.plot(
                    result.times,
                    result.expect[i],
                    linewidth=1,
                    label="{}".format(i + 1),
                )
        ax.legend()
        ax.set_ylim([-0.01, 1.01])
        ax.set_xlabel("Time")
        ax.set_ylabel("Population")
        plt.show()

    def _plot_population_diagonalization(self, tlist, result, dim):
        fig, ax = plt.subplots()
        tlist = np.array(tlist)
        n = dim
        for i in range(dim):
            ax.plot(
                tlist, np.real(result[:, i]), linewidth=1.5, label="{}".format(i + 1)
            )
        ax.legend()
        ax.set_ylim([-0.1, 1.1])
        ax.set_xlabel(r"Time")
        ax.set_ylabel(r"Population")
        plt.show()

    def simulate(
        self,
        initial_state_index_list,
        level_to_print_max_value,
        maxtime,
        delta_stark_shift=0,
        Diagonalization=True,
        plot_pop=True,
        Trotterintervals=500,
        points_per_TI=2,
        resolution=250,
    ):
        """
        A function to simulate the time evolution of the population of every level. It also returns the maximum population of specifically the level_to_print_max_value.

        Args:
            initial_state_index_list (list): value for :attr:`initial_state_index_list`
            level_to_print_max_value (number): value for :attr:`level_to_print_max_value`
            maxtime (integer): value for :attr:`maxtime`
            delta_stark_shift (number): value for :attr:`delta_stark_shift`
            Diagonalization (bool): value for :attr:`Diagonalization`
            plot_pop (bool): value for :attr:`plot_pop`
            Trotterintervals (number): value for :attr:`Trotterintervals`
            points_per_TI (number): value for :attr:`points_per_TI`
            resolution (number): value for :attr:`resolution`

        Attributes:
            initial_state_index_list (list): Initial population distribution. First entry denotes the population of the first level and so on. Length of the list needs to be equal to the level-count and the entries need to sum up to one.
            level_to_print_max_value (number): Just affects the output. Determines which levels population maximum is printed in the output. 0 is the first level.
            maxtime (integer): Maximum time the system is simulated.
            delta_stark_shift (number): Detuning of the first level.
            Diagonalization (bool): If False, simulate the system with the integration method 'qutip.mesolve' from QuTiP. If True, simulate the system using diagonalization.
            plot_pop (bool): If False, do not show the population plot.
            Trotterintervals (number): Only relevant if a pulse is given. Discretizes the pulse into a step function of `Trotterintervals` time intervals.
            points_per_TI (number): Only relevant if a pulse is given. Divides one trotterinterval in points_per_TI intervals. This determines the number of points calculated within one trotterinterval.
            resolution (number): Only relevant if Diagonalization = True and no pulse is given. Divides the simulation time interval into `resolution` uniformly distributed points of time.

        """
        Trotter = False
        # make basis: levels to kets
        initial_state_dm = qutip.qzero(self.dim)
        level_ket = [[]] * self.dim  # empty ket with 'dim' empty entries
        proj = [[]] * self.dim
        for i in range(len(self.levels)):
            level_ket[i] = qutip.basis(self.dim, i)
            proj[i] = qutip.ket2dm(level_ket[i])
            initial_state_dm = (
                initial_state_dm
                + qutip.ket2dm(level_ket[i]) * initial_state_index_list[i]
            )
        # level_ket: [Level, corresponding basis vector], e.g. [[level1, level2], [[1,0],[0,1]]]
        level_ket = [
            self.levels,
            level_ket,
        ]
        # generation of the ladder operators:
        sigma = [[]] * len(
            self.decay.rates
        )  # for each decay rate there is one ladder operator
        start_state = [[]] * len(
            self.decay.rates
        )  # for each gamma there is a starting state and a target state
        final_state = [[]] * len(self.decay.rates)
        for i in range(len(self.decay.rates)):
            if (
                self.decay.final_states[i][0].energy
                > self.decay.final_states[i][1].energy
            ):
                start_state = self.decay.final_states[i][0]
                final_state = self.decay.final_states[i][1]
            else:
                start_state = self.decay.final_states[i][1]
                final_state = self.decay.final_states[i][0]
            i1 = level_ket[0].index(start_state)
            i2 = level_ket[0].index(final_state)
            sigma[i] = (
                level_ket[1][i2] * level_ket[1][i1].dag()
            )  # sigma[i] is the ladder operator for the first decay rate in the list of the decay object
        rabifreqs = [[]] * len(self.lasers)

        # Construction of the Hamiltonian H ---------
        # DIAGONAL ENTRIES:
        E_virt = [[]] * len(self.levels)
        Bool_Detuning = [[]] * len(self.levels)
        i1 = [[]] * len(self.lasers)
        i2 = [[]] * len(self.lasers)
        for i in range(len(self.levels)):
            if i == 0:
                E_virt[i] = 0 + delta_stark_shift
            else:
                for n in range(len(self.lasers)):
                    i1[n] = level_ket[0].index(
                        self.lasers[n].couple[0]
                    )  # level index of the level in couple[0]. (e.g. index of level1 in couple [level1, level3])
                    i2[n] = level_ket[0].index(
                        self.lasers[n].couple[1]
                    )  # level index of the level in couple[1]. (e.g. index of level3 in couple [level1, level3])
                for n in range(len(self.lasers)):
                    if (
                        self.lasers[n].couple[1] == self.levels[i]
                    ):  # If the target entry in the lasercouple is equal to the level with E_virt[i].
                        if Bool_Detuning[i] == True:
                            if Bool_Detuning[i1[n]] == True:
                                if np.around(
                                    E_virt[i] - E_virt[i1[n]], decimals=3
                                ) == np.around(self.lasers[n].frequency, decimals=3):
                                    pass
                                else:
                                    print(E_virt[i] - E_virt[i1[n]])
                                    print(self.lasers[n].frequency)
                                    print(
                                        np.around(E_virt[i] - E_virt[i1[n]], decimals=4)
                                    )
                                    print(
                                        np.around(self.lasers[n].frequency, decimals=4)
                                    )
                                    raise TypeError(
                                        "Not enough degrees of freedom: See my Master Thesis page 7 bottom text: The only prerequisite ..."
                                    )
                            else:
                                E_virt[i1[n]] = E_virt[i2[n]] - self.lasers[n].frequency
                                Bool_Detuning[i1[n]] == True
                        else:
                            E_virt[i] = E_virt[i1[n]] + self.lasers[n].frequency
                            Bool_Detuning[
                                i
                            ] = True  # This level does not support further detuning and further detunings must be "passed on" to other levels.
                    else:
                        if Bool_Detuning[i] == True:
                            pass
                        else:
                            E_virt[i] = self.levels[
                                i
                            ].energy  # No detuning. This can change in the following steps, if you have to "pass on" a detuning to this level.

        H = np.zeros((self.dim, self.dim), dtype=np.complex128)
        for i in range(self.dim):
            for j in range(self.dim):
                if i == j:
                    H[i, j] = E_virt[i] - self.levels[i].energy  # = detuning
        # NON-DIAGONAL ENTRIES:
        H_couple1 = [[]] * len(self.lasers)
        H_couple2 = [[]] * len(self.lasers)
        for i in range(len(self.lasers)):
            rabifreqs[i] = self.lasers[i].rabi
            i1 = level_ket[0].index(
                self.lasers[i].couple[0]
            )  # level index of the level in couple[0]. (e.g. index of level1 in couple [level1, level3])
            i2 = level_ket[0].index(
                self.lasers[i].couple[1]
            )  # level index of the level in couple[1]. (e.g. index of level3 in couple [level1, level3])
            basis1 = level_ket[1][
                i1
            ]  # gives the corresponding basis vector to i1, i.e. to couple[0].
            basis2 = level_ket[1][i2]  # gives the corresponding basis vector to i2
            H_couple1[i] = (
                basis1 * basis2.dag()
            )  # gives the right place in the H matrix for the first entry of the coupling Rabi frequency of the i-th laser
            H_couple2[i] = (
                basis2 * basis1.dag()
            )  # gives the right place in the H matrix for the second entry. H_couple1[i] and H_couple2[i] form a pair for the i-th laser

        # Time-dependent pulses
        lasermitpuls_index = []
        for n in range(
            len(self.lasers)
        ):  # Here it is decided whether Trotter should be performed (i.e. it is checked whether a pulse is given or not)
            if hasattr(
                self.lasers[n].pulse, "__len__"
            ):  # If the pulse is given as a list
                if self.lasers[n].pulse.any() != None:  # If laser[n] has a pulse
                    lasermitpuls_index = np.append(
                        lasermitpuls_index, n
                    )  # All laser indices from the lasers that have a pulse are in here
                    lasermitpuls_index = lasermitpuls_index.astype(int)
                    Trotter = True
                else:
                    H = H + (1 / 2 * rabifreqs[n]) * (H_couple1[n] + H_couple2[n])
            else:  # If the pulse is given as a function
                if self.lasers[n].pulse != None:  # If laser[n] has a pulse
                    lasermitpuls_index = np.append(
                        lasermitpuls_index, n
                    )  # All laser indices from the lasers that have a pulse are in here
                    lasermitpuls_index = lasermitpuls_index.astype(int)
                    Trotter = True
                else:
                    H = H + (1 / 2 * rabifreqs[n]) * (H_couple1[n] + H_couple2[n])
        print(
            "Hamiltonian in the rotating frame: {}".format(H)
        )  # Only the Rabi frequency entries are missing if the pulse is time dependent.

        # The Hamiltonian is created, now we solve it.
        # collapse operators are defined
        c_ops = [[]] * len(self.decay.rates)
        for i in range(len(self.decay.rates)):
            c_ops[i] = qutip.Qobj(np.sqrt(self.decay.rates[i]) * sigma[i])
        # !Entweder Diagonalisieren oder mesolve verwenden!
        if (
            Trotter == True
        ):  # Trotter intervals: Intervals in which H is assumed to be constant. Time intervals: Intervals for which a value is calculated (important for the balance between resolution and calculation time, determines how many values are calculated in total).
            number_TI = Trotterintervals  # number of trotterintervals
            if maxtime / number_TI != int(maxtime / number_TI):
                raise Exception(
                    "The number of trotterintervals of {} doesnt divide the given timerange of {} into trotterintervals with integer time".format(
                        number_TI, maxtime
                    )
                )
            td = np.linspace(0, maxtime, number_TI + 1)
            print("One trotterinterval has size {}.".format(td[1]))
            if td[1] / points_per_TI != int(td[1] / points_per_TI):
                raise Exception(
                    "The number of timeintervals per trotterinterval of {} doesnt divide one trotterinterval of size {} into timeintervals with integer time".format(
                        points_per_TI, td[1]
                    )
                )
            print("One trotter step has size {}.".format(td[1] / points_per_TI))
            trotter_step = int(td[1] / points_per_TI)

            # Time evolution
            tlist_ges = [[]] * (
                (len(td) - 1) * int((maxtime / (number_TI)) / trotter_step) + 1
            )  # will contain every point of time
            dm_t = [[]] * len(tlist_ges)  # will contain density matrices for each time
            result = [[]] * len(
                tlist_ges
            )  # will contain one list of expectation values for each level for each time

            # For t = 0, initial population:
            expect_value_0 = np.zeros(self.dim)
            tlist_ges[0] = 0
            dm_t[0] = np.reshape(initial_state_dm, (self.dim**2, 1))
            for m in range(self.dim):
                expect_value_0[m] = qutip.expect(
                    qutip.Qobj(proj[m]),
                    qutip.Qobj(np.reshape(dm_t[0], (self.dim, self.dim))),
                )
            result[0] = np.array(expect_value_0)

            H_list = [H] * len(td)
            L = [[]] * len(td)
            # For t > 0:
            for n in range(len(td) - 1):  # loop over all trotter intervals
                # add the Rabi frequency entries for the Trotter interval we want to simulate here
                for i in range(len(lasermitpuls_index)):
                    if hasattr(
                        self.lasers[0].pulse, "__len__"
                    ):  # If the pulse is given as a list
                        H_list[n] = H_list[n] + (
                            1 / 2 * self.lasers[lasermitpuls_index[i]].pulse[n]
                        ) * (
                            H_couple1[lasermitpuls_index[i]]
                            + H_couple2[lasermitpuls_index[i]]
                        )
                    else:  # If the pulse is given as a function
                        H_list[n] = H_list[n] + (
                            1 / 2 * self.lasers[lasermitpuls_index[i]].pulse(td[n])
                        ) * (
                            H_couple1[lasermitpuls_index[i]]
                            + H_couple2[lasermitpuls_index[i]]
                        )
                L[n] = qutip.liouvillian(H_list[n], c_ops)
                L[n] = np.reshape(L[n], (self.dim**2, self.dim**2))
                # Time evolution within one trotter interval:
                tlist = list(
                    range(
                        int(td[n]),
                        int(td[n + 1] + 1),
                        trotter_step,
                    )
                )
                t_index = range(0, len(tlist))
                # For the first trotter interval:
                if (
                    n == 0
                ):  # +n*len(tlist) braucht man immer, wenn es eine "globale" Variable ist und keine, die sich in jedem Zykel selbst überschreiben soll.   # tlist[1] = t_index[1]+n*len(tlist)
                    for t in range(t_index[1], t_index[-1] + 1):
                        tlist_ges[t] = tlist[t]
                        expect_value = [
                            []
                        ] * self.dim  # will contain expectation values of all levels for one t
                        dm_t[t] = scipy.linalg.expm(L[n] * tlist[t]) @ np.reshape(
                            dm_t[0], (self.dim**2, 1)
                        )
                        for m in range(self.dim):
                            expect_value[m] = qutip.expect(
                                qutip.Qobj(proj[m]),
                                qutip.Qobj(np.reshape(dm_t[t], (self.dim, self.dim))),
                            )
                        result[t] = np.array(expect_value)
                # For all other trotter intervals:
                else:
                    for t in range(t_index[1], t_index[-1] + 1):
                        tlist_ges[t + n * points_per_TI] = tlist[t]
                        expect_value = [
                            []
                        ] * self.dim  # will contain expectation values of all levels for one t
                        dm_t[t + n * points_per_TI] = scipy.linalg.expm(
                            L[n] * (tlist[t] - tlist[0])
                        ) @ np.reshape(dm_t[n * points_per_TI], (self.dim**2, 1))
                        for m in range(self.dim):
                            expect_value[m] = qutip.expect(
                                qutip.Qobj(proj[m]),
                                qutip.Qobj(
                                    np.reshape(
                                        dm_t[t + n * points_per_TI],
                                        (self.dim, self.dim),
                                    )
                                ),
                            )
                        result[t + n * points_per_TI] = np.array(expect_value)
            result = np.array(
                result
            )  # result[t] is a list of expectation values of all levels for the point of time that is indexed with t
            if plot_pop == True:
                self._plot_population_diagonalization(tlist_ges[:], result[:], self.dim)
            print(
                "Maximum population of level {}:".format(level_to_print_max_value + 1)
            )
            return np.real(np.amax(result[:, level_to_print_max_value]))

        elif Diagonalization == True:
            # Liouvillian
            L = qutip.liouvillian(H, c_ops)
            L = np.reshape(L, (self.dim**2, self.dim**2))

            # Time evolution:
            tlist = np.linspace(
                0, maxtime, resolution
            )  # the time interval is divided into 250 points of time
            dm_t = [[]] * len(tlist)  # will contain density matrices for each time
            result = [[]] * len(
                tlist
            )  # will contain one list of expectation values for each level for each time

            # For t = 0:
            expect_value_0 = [
                []
            ] * self.dim  # will contain expectation values for every level for t = 0
            dm_t[0] = np.reshape(initial_state_dm, (self.dim**2, 1))
            for m in range(self.dim):
                expect_value_0[m] = qutip.expect(
                    qutip.Qobj(proj[m]),
                    qutip.Qobj(np.reshape(dm_t[0], (self.dim, self.dim))),
                )
            result[0] = np.array(expect_value_0)

            # For t > 0:
            for t in range(1, len(tlist)):
                expect_value = [
                    []
                ] * self.dim  # will contain expectation values of all levels for one t
                dm_t[t] = scipy.linalg.expm(L * tlist[t]) @ np.reshape(
                    dm_t[0], (self.dim**2, 1)
                )
                for m in range(self.dim):
                    expect_value[m] = qutip.expect(
                        qutip.Qobj(proj[m]),
                        qutip.Qobj(np.reshape(dm_t[t], (self.dim, self.dim))),
                    )
                result[t] = np.array(np.reshape(expect_value, (self.dim,)))
            result = np.array(
                result
            )  # result[t] is a list of expectation values of all levels for the point of time that is indexed with t
            if plot_pop == True:
                self._plot_population_diagonalization(tlist, result, self.dim)
            print(
                "Maximum population of level {}:".format(level_to_print_max_value + 1)
            )
            return np.real(np.amax(result[:, level_to_print_max_value]))

        else:  # Integration-solver from QuTip. See QuTiP Documentation of the mesolve function.
            tlist = list(range(0, maxtime, int(41300000 / 100)))
            print("Length of tlist: {}".format(len(tlist)))
            # opts=Options(nsteps=1000)
            # print(options)
            # start = time.time()
            result = qutip.mesolve(H, initial_state_dm, tlist, c_ops, proj)
            # end = time.time()
            # print(f"{end-start:5.3f}s")
            self._plot_population(result, self.dim)
            print(result.expect[1][-1])
