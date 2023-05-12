#%%
# Imports
import numpy as np
import qutip
import matplotlib
import matplotlib.pylab as plt
import scipy
import time

def plot_population(result, dim, maxdetuning=None):
    fig, ax = plt.subplots()
    for i in range(dim):
        if i != 12:
            ax.plot(result.times, result.expect[i],linewidth = 1, label='{}'.format(i+1))
            #ax.plot(result.times*2.41888e-17 , result.expect[i],linewidth = 1, label='{}'.format(i+1))
#    ax.plot(result.times, result.expect[0]+result.expect[1]+result.expect[2], label='total')
    # x = np.arange(0,50000,10)
    # y = np.sin(1.2441522993618834e-05 /2*x)**2     # sin(Omega*x) sollte es ein: "The Rabi frequency is the radian frequency of the Rabi cycle undergone for a given atomic transition in a given light field" <==> sin**2(Omega/2*x)
    # ax.plot(x,y, linewidth=3)
    # x2 = np.arange(0,200000,100)
    # y2 = np.sin(4.9766091974475336e-05/2*x2)**2
    # ax.plot(x2,y2, linewidth=1)
    # if maxdetuning != None:
    #     ax.vlines((maxdetuning), 0,1)
    ax.legend()
    ax.set_ylim([-0.01, 1.01])
    ax.set_xlabel('Time')
    # ax.set_xlabel('time in s')
    ax.set_ylabel('Population')
    # plt.savefig('bigplot.png')
    plt.show()

def plot_population_diagonalization(tlist, result, dim):
    fig, ax = plt.subplots()
    tlist = np.array(tlist)
    n = dim
    for i in range(dim):
        ax.plot(tlist/41341373000, np.real(result[:,i]) ,linewidth = 1.5)
    ax.legend()
    ax.set_ylim([-0.1, 1.1])
    ax.set_xlabel(r'Time in µs')
    ax.set_ylabel(r'Population')
    plt.show()

class Level:
    '''
    Level characterized by the energy
    '''
    def __init__(self, energy):
        if type(energy) != list:
            raise TypeError('energy needs to be a list') # ist inzwischen unnötig das als Liste zu machen
        self.energy = energy[0]

class Decay:
    '''Ein Objekt gibt alle decays an. Fuer jedes couple in final_states muss es ein gamma in rates geben.'''
    def __init__(self, rates, final_states): # final_states is a Level-list with two-level-couples in each entry. rates is the respective gamma-list to each couple
        if (type(rates) != list or type(final_states) != list):
            raise TypeError('rates and final_states need to be a list')
        self.rates = np.array(rates)
        self.final_states = np.array(final_states)

class Laser:
    '''
    characterizing values for the laser: Rabi-Frequency, frequency and coupled states.
    '''
    def __init__(self, rabifreq, frequency, couple, polarization=None, pulse=None): # couple is a list of two Levels e.g. [level1, level2]
        self.rabi = rabifreq  # Omega komplexwertig ist die Rabifrequenz
        self.couple = couple
        self.frequency = frequency  # omega
        self.polarization = polarization # Das ist eine Liste aus einem normierten E-Feld Vektor im Laserkoordinatensystem, einem Theta_k und einem Theta_p (in Grad) (siehe Skizze auf dem Blatt "Polarization") und einem Paar [m_i,m_f]
        self.pulse = pulse # Das ist die Puls-Funktion, also die zeitabhängige Rabi-Frequenz. In der Form einer Funktion "pulse(t)"

class System:
    '''
    methods: draw and get_population
    lasers & levels: list of lasers & levels that act in the system
    '''
    def __init__(self, levels, lasers, decay):  # levels soll Liste von Level-Objekten sein. !Reihenfolge bestimmt Basisvektoren!
        self.lasers = lasers
        self.levels = levels
        self.dim = len(levels)   # dim gibt die Anzahl an Leveln
        self.decay = decay

    def draw(self):
        #x_window = np.linspace(0,1,20)     # self.levels -> System.levels -> levels (das ist eine Liste von Level-Objekten)
        # Levels
        for i in range(len(self.levels)):   # draw hat keinen Parameter levels, aber das System-Objekt hat die Methode .levels aus __init__
            x = [0.01+0.2*i,0.19+0.2*i]
            y = [self.levels[i].energy,self.levels[i].energy]
            plt.plot(x, y, linestyle='-',linewidth=1, label='Level {}'.format(i+1))  # jedes Level-Objekt hat die Methode .energy
            plt.text(0.5*(x[0]+x[1]), y[0]+0.002, '{}'.format(i+1))
        plt.xlim([-0.25,1.5])
        # Detuning. Lasers tuned to a frequency below the resonant frequency are called 'red-detuned' andere Seite ist 'blue-detuned'.
        detuning = [[]] * len(self.lasers)
        for i in range(len(self.lasers)):
            max_couple = max(self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy)
            min_couple = min(self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy)
            detuning[i] = min_couple + self.lasers[i].frequency - max_couple # freq zu klein -> negatives detuning. hbar=1
            if abs(detuning[i]) > 2: # self.lasers[i].couple[0].energy + self.lasers[i].frequency - self.lasers[i].couple[1].energy < self.lasers[i].couple[0].energy - self.lasers[i].frequency - self.lasers[i].couple[1].energy
                detuning[i] = - (self.lasers[i].couple[0].energy - self.lasers[i].frequency - self.lasers[i].couple[1].energy)
        for i in range(len(self.lasers)):
            max_couple = max(self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy)
            if max_couple == self.lasers[i].couple[0].energy:
                i1 = self.levels.index(self.lasers[i].couple[0])
            else:
                i1 = self.levels.index(self.lasers[i].couple[1])
            x = [0.01+0.2*i1,0.19+0.2*i1]
            y = [max_couple + detuning[i],max_couple + detuning[i]]
            plt.plot(x, y, linestyle='--', label='D {} = {:.2e}'.format(i+1, detuning[i]/Hz_to_au))
        # Lasers
        for i in range(len(self.lasers)):
            max_couple = max(self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy)
            min_couple = min(self.lasers[i].couple[0].energy, self.lasers[i].couple[1].energy)
            i1 = self.levels.index(self.lasers[i].couple[0])
            i2 = self.levels.index(self.lasers[i].couple[1])
            x = [0.1+0.2*i1,0.1+0.2*i2]
            # x = [0.1+i/len(self.lasers),0.1+i/len(self.lasers)]
            y = [min_couple,max_couple + detuning[i]]
            plt.plot(x, y, linestyle='-', label='Laser {}'.format(i+1))
            # plt.text(x[0]+0.01, 0.5*(y[0]+y[1]), '{}'.format(i+1))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()


    def fidelity(self, initial_state_index_list, read_out_level, maxtime, delta_stark_shift, Diagonalization=True, plot_pop=True):  # Level nach Energie aufsteigend sortieren. Laserkoppelpaare von niedrig nach groß.
        '''
        hi
        '''
        Trotter = False
        # make basis: levels to kets
        initial_state_dm = qutip.qzero(self.dim)
        # initial_state = basis(self.dim, initial_state_index)
        level_ket = [[]] * self.dim  # empty ket with dim empty entries
        proj = [[]] * self.dim
        for i in range(len(self.levels)):
            level_ket[i] = qutip.basis(self.dim, i)
            proj[i] = qutip.ket2dm(level_ket[i])
            # print(proj[i])
            initial_state_dm = initial_state_dm + qutip.ket2dm(level_ket[i])*initial_state_index_list[i]
        level_ket = [self.levels,level_ket]     # Paare: [Level, zugehoeriger Basisvektor] e.g. [[level1, level2], [[1,0],[0,1]]]
        # Erzeugung der Leiteroperatoren:
        sigma = [[]] * len(self.decay.rates) # fuer jedes Gamma gibt es einen Leiteroperator
        start_state = [[]] * len(self.decay.rates) # fuer jedes Gamma gibt es einen Start-
        final_state = [[]] * len(self.decay.rates) # und einen Zielzustand
        for i in range(len(self.decay.rates)):
            if (self.decay.final_states[i][0].energy > self.decay.final_states[i][1].energy):
                start_state = self.decay.final_states[i][0]
                final_state = self.decay.final_states[i][1]
            else:
                start_state = self.decay.final_states[i][1]
                final_state = self.decay.final_states[i][0]
            i1 = level_ket[0].index(start_state)
            i2 = level_ket[0].index(final_state)
            sigma[i] = level_ket[1][i2] * level_ket[1][i1].dag()    # sigma[i] gibt Leiteroperator fuer das erste Gamma in der Liste vom decay-Objekt
        rabifreqs = [[]] * len(self.lasers)
        freqs = [[]] * len(self.lasers)
        ############## U Calculations
        U = np.zeros((self.dim, self.dim), dtype=np.complex128)
        U_exp = [[]] * self.dim

        # DIAGONALEINTRAEGE:
        # detuning = [[]] * len(self.levels)
        E_virt = [[]] * len(self.levels)
        Bool_Detuning = [[]] * len(self.levels)
        # Vergeben = [[]] * len(self.levels)
        i1 = [[]] * len(self.lasers)
        i2 = [[]] * len(self.lasers)
        for i in range(len(self.levels)):
            # print(i)
            if i == 0:
                E_virt[i] = 0 + delta_stark_shift   # Man koennte evtl noch einen Freiheitsgrad rausschlagen, wenn man hier auch ein Detuning zulassen wuerde. Aber andererseits ist das dann aequivalent zu einem Energieshift aller anderen Level, die in diesem Algorithmus dann auch auftauchen sollten, also ist es vermutlich doch nicht moeglich.
            else:
                for n in range(len(self.lasers)):
                    i1[n] = level_ket[0].index(self.lasers[n].couple[0]) # gibt Level-Index vom Level in couple[0] (Index von level1 bei couple [level1, level3])
                    i2[n] = level_ket[0].index(self.lasers[n].couple[1]) # gibt Level-Index vom Level in couple[1] (Index von level3 bei couple [level1, level3])
                # print(i1)
                for n in range(len(self.lasers)):
                    # print(n)
                    if self.lasers[n].couple[1] == self.levels[i]:  # Wenn der Zieleintrag im Lasercouple gleich dem Level mit E_virt[i] ist. ########
                        if Bool_Detuning[i] == True:
                            if Bool_Detuning[i1[n]] == True:
                                if np.around(E_virt[i] - E_virt[i1[n]], decimals=3) == np.around(self.lasers[n].frequency, decimals=3):
                                    pass
                                else:
                                    print(E_virt[i] - E_virt[i1[n]])
                                    print(self.lasers[n].frequency)
                                    print(np.around(E_virt[i] - E_virt[i1[n]], decimals=4))
                                    print(np.around(self.lasers[n].frequency, decimals=4))
                                    raise TypeError('Not enough degrees of freedom: Laserfrequenzen sorgen dafuer, dass ein Level mehrere Energien annehmen muss bzw. mindestens zwei verschiedene Detunings besitzt. See my Master Thesis page 7 bottom text: the only prerequisite ...')
                            else:
                                E_virt[i1[n]] = E_virt[i2[n]] - self.lasers[n].frequency ##########
                                Bool_Detuning[i1[n]] == True
                        else:
                            E_virt[i] = E_virt[i1[n]] + self.lasers[n].frequency   # E2v = E1v + omega_12
                            # Vergeben[i] == True
                            Bool_Detuning[i] = True  # Dieses Level vertraegt kein weiteres Detuning und weitere Detunings muessen an andere Level "weitergegeben" werden
                    else: #############
                        if Bool_Detuning[i] == True:
                            pass
                        else:
                            E_virt[i] = self.levels[i].energy   # Kein Detuning. Kann sich aber noch in folgenden Schritten aendern, wenn man auf dieses Level ein Detuning "weitergeben" muss. siehe Beispiel in "Neuer ALGO Seite 4".

        H = np.zeros((self.dim,self.dim), dtype=np.complex128)
        for i in range(self.dim):
            for j in range(self.dim):
                if i == j:
                    H[i,j] = E_virt[i] - self.levels[i].energy # = detuning
        # NICHT-DIAGONAL EINTRAEGE:
        H_couple1 = [[]] * len(self.lasers)
        H_couple2 = [[]] * len(self.lasers)
        for i in range(len(self.lasers)):
            rabifreqs[i] = self.lasers[i].rabi
            i1 = level_ket[0].index(self.lasers[i].couple[0])  # gibt Level-Index vom Level in couple[0] (Index von level1 bei couple [level1, level3])
            i2 = level_ket[0].index(self.lasers[i].couple[1])  # gibt Level-Index vom Level in couple[1] (Index von level3 bei couple [level1, level3])
            basis1 = level_ket[1][i1]   # gibt zugehoerigen Basisvektor zu i1, also zu couple[0]
            basis2 = level_ket[1][i2]   # gibt zugehoerigen Basisvektor zu i2
            H_couple1[i] = basis1 * basis2.dag() # gibt den Platz fuer den ersten Eintrag des Kopplungs-Omegas des i-ten Lasers
            H_couple2[i] = basis2 * basis1.dag() # gibt den Platz fuer den zweiten Eintrag. H_couple1[i] und H_couple2[i] bilden ein Paar fuer i-ten Laser

        # zeitabhängige Pulse rein
        lasermitpuls_index = []
        for n in range(len(self.lasers)): # Hier wird entschieden, ob Trotter gemacht werden soll (d.h. es wird geschaut, ob ein Puls gegeben ist oder nicht)
            if hasattr(self.lasers[n].pulse, "__len__"): # Falls der Puls als Liste gegeben ist
                if self.lasers[n].pulse.any() != None: # Falls laser[n] einen Puls gegeben hat
                    lasermitpuls_index = np.append(lasermitpuls_index, n) # Alle Laserindizes von den Lasern, die einen Puls haben, sind hier drin
                    lasermitpuls_index = lasermitpuls_index.astype(int)
                    Trotter = True
                    # print(lasermitpuls_index)
                else:
                    H =  H + (1/2*rabifreqs[n]) * (H_couple1[n] + H_couple2[n])
            else: # Falls der Puls als Funktion gegeben ist
                if self.lasers[n].pulse != None: # Falls laser[n] einen Puls gegeben hat
                    lasermitpuls_index = np.append(lasermitpuls_index, n) # Alle Laserindizes von den Lasern, die einen Puls haben, sind hier drin
                    lasermitpuls_index = lasermitpuls_index.astype(int)
                    Trotter = True
                    # print(lasermitpuls_index)
                else:
                    H =  H + (1/2*rabifreqs[n]) * (H_couple1[n] + H_couple2[n])

        # Hier an dieser Stelle fehlen nur noch die Rabifrequenzen der Pulse.

        # print('Alles ist transformiert: {}' .format(H))
        # Loesen
        c_ops = [[]] * len(self.decay.rates)
        for i in range(len(self.decay.rates)):
            c_ops[i] = qutip.Qobj(np.sqrt(self.decay.rates[i])*sigma[i])
        # !Entweder Diagonalisieren oder mesolve verwenden!
        if Trotter == True: # Trotterintervalle: Intervalle, in denen H konstant angenommen wird. Zeitintervalle: Intervalle für die ein Wert ausgerechnet wird (Wichtig für die Balance zwischen Auflösung und Rechenzeit, bestimmt wie viele Werte insgesamt berechnet werden).
            Anzahl_n = 500+1 # Anzahl_n-1 = Anzahl Trotterintervalle. Bestimmt die Größe der Trotterintervalle, auf denen der Hamiltonian/Liouvillian konstant angenommen wird. Anzahl_n = len(zeitdiskretisierung). Je größer, desto mehr Trotterintervalle. 
            if maxtime/(Anzahl_n-1) != int(maxtime/(Anzahl_n-1)):
                raise Exception('The number of trotterintervals of {} doesnt divide the given timerange of {} into trotterintervals with integer time'.format(Anzahl_n-1, maxtime))
            zeitdiskretisierung = np.linspace(0,maxtime,Anzahl_n) # macht "Anzahl_n" Einträge von 0 bis maxtime. !!!!!!!!!!NUR GANZZAHLIGE WERTE IN DIESER LISTE ERLAUBT!!!!!!!!!!, sonst könnte das Runden unten bei tlist zu Problemen mit den Dimensionen führen
            print('Ein Trotterintervall ist {} groß.'.format(zeitdiskretisierung[1]))
            anzahl_zeitintervalle_pro_trotterintervall = 2 # muss größer als 1 sein und ein Trotterintervall in integer Zeitintervalle teilen.
            if zeitdiskretisierung[1]/anzahl_zeitintervalle_pro_trotterintervall != int(zeitdiskretisierung[1]/anzahl_zeitintervalle_pro_trotterintervall):
                raise Exception('The number of timeintervals per trotterinterval of {} doesnt divide one trotterinterval of size {} into timeintervals with integer time'.format(anzahl_zeitintervalle_pro_trotterintervall, zeitdiskretisierung[1]))
            print(zeitdiskretisierung[1]/anzahl_zeitintervalle_pro_trotterintervall)
            trotter_step = int(zeitdiskretisierung[1]/anzahl_zeitintervalle_pro_trotterintervall)   #41341373*20 # Muss ganzzahlig sein. # hat nichts mit Genauigkeit zu tun, sondern nur mit der "Auflösung" des Bildes bzw. Anzahl ausgerechneter Werte. Ein Trotter_step gibt an wie groß das Intervall zwischen ausgerechneten Werten ist: Je größer, desto weniger Werte werden ausgerechnet.
            # Der Code ist schlecht, deshalb muss trotter_step erstens kleiner sein als ein Trotterintervall und zweitens so gewählt werden, dass die Trotterintervalle ganzzahlig zerlegt werden. D.h. trotter_step muss zeitdiskretisierung[1] ganzzahlig teilen.
            # Für t=0 bei der Diagonalisierung:
            tlist_ges = [[]] * ((len(zeitdiskretisierung)-1) * int((maxtime/(Anzahl_n-1))/trotter_step)+1)    # len(zeitdiskretisierung-1) ist die Anzahl der Trotter-Intervalle.  maxtime/(Anzahl_n-1) ist die Größe eines Trotterintervalls in Zeiteinheiten.   (maxtime/(Anzahl_n-1))/trotter_step ist die Anzahl an Zeitpunkten, die berechnet werden soll innerhalb eines Trotterintervalls (also innerhalb einem n)
            print((len(zeitdiskretisierung)-1) * int((maxtime/(Anzahl_n-1))/trotter_step))
            dm_t = [[]] * len(tlist_ges) # Eine Liste, die für jeden Zeitpunkt die zugehörige Dichtematrix enthalten wird
            result = [[]] * len(tlist_ges)
            zwischenresult0 = np.zeros(self.dim, dtype=np.complex128) #[[]] * self.dim
            tlist_ges[0] = 0
            dm_t[0] = np.reshape(initial_state_dm, (self.dim**2,1))
            # dm_t[0] = np.reshape(proj[initial_state_index], (self.dim**2,1))
            for m in range(self.dim): # für t=0 werden die expectation values zu jedem Level ausgerechnet
                zwischenresult0[m] = qutip.expect(qutip.Qobj(proj[m]), qutip.Qobj(np.reshape(dm_t[0], (self.dim,self.dim))))
            result[0] = np.array(zwischenresult0, dtype=np.complex128)
            # print(self.lasers[0].pulse(1))
            print(len(zeitdiskretisierung))
            
            H_list = [H] * len(zeitdiskretisierung)
            L = [[]] * len(zeitdiskretisierung)
            # pulse = np.array([pulse(t) for t in zeitdiskretisierung]) # Gibt den Wert der Rabifrequenz pro Zeitintervall im Verlauf des Pulses als Liste
            for n in range(len(zeitdiskretisierung)-1): # n loopt über die Intervalle
                if n % 100 == 0:
                    print(n)
                for i in range(len(lasermitpuls_index)):
                    if hasattr(self.lasers[0].pulse, "__len__"): # Falls der Puls als Liste gegeben ist
                        H_list[n] = H_list[n] + (1/2*self.lasers[lasermitpuls_index[i]].pulse[n]) * (H_couple1[lasermitpuls_index[i]] + H_couple2[lasermitpuls_index[i]])
                    else: # Falls der Puls als Funktion gegeben ist
                        H_list[n] = H_list[n] + (1/2*self.lasers[lasermitpuls_index[i]].pulse(zeitdiskretisierung[n])) * (H_couple1[lasermitpuls_index[i]] + H_couple2[lasermitpuls_index[i]])
                L[n] = qutip.liouvillian(H_list[n], c_ops)
                L[n] = np.reshape(L[n], (self.dim**2, self.dim**2))
                # Vollständige Zeitentwicklung:
                tlist = list(range(int(zeitdiskretisierung[n]),int(zeitdiskretisierung[n+1]+1), trotter_step)) # Für den Zeitraum, in dem ein L[n] gilt.  int(zeitdiskretisierung[n+1]) = int(zeitdiskretisierung[n]) + maxtime/(Anzahl_Steps-1)
                t_index = range(0, len(tlist)) # ist identisch zu tlist für n = 0, wenn trotter_step = 1
                if n == 0:  #+n*len(tlist) braucht man immer, wenn es eine "globale" Variable ist und keine, die sich in jedem Zykel selbst überschreiben soll.   # tlist[1] = t_index[1]+n*len(tlist)
                    for t in range(t_index[1],t_index[-1]+1): # Loopt über die Anzahl an Zeitpunkten innerhalb einer konstanter-Puls-Phase. Jeder Zeitpunkt wird mit t geINDEXT
                        tlist_ges[t] = tlist[t]
                        zwischenresult = [[]] * self.dim # expectation values für ein t zu jedem Level
                        dm_t[t] = scipy.linalg.expm(L[n]*tlist[t]) @ np.reshape(dm_t[0], (self.dim**2,1)) 
                        for m in range(self.dim): # für jedes t werden die expectation values zu jedem Level ausgerechnet
                            # if m == 1:
                            #     print(Qobj(np.reshape(dm_t[t], (self.dim,self.dim))))
                            zwischenresult[m] = qutip.expect(qutip.Qobj(proj[m]), qutip.Qobj(np.reshape(dm_t[t], (self.dim,self.dim))))
                        result[t] = np.array(zwischenresult, dtype=np.complex128)   # für jedes t werden die 7 (je nach levelanzahl) expectation values hier eingesetzt. D.h. result[t] gibt eine Liste mit den Expectation values aller Level
                else:
                    for t in range(t_index[1],t_index[-1]+1): # Loopt über die Anzahl an Zeitpunkten. Jeder Zeitpunkt wird mit t geindext
                        tlist_ges[t+n*anzahl_zeitintervalle_pro_trotterintervall] = tlist[t]  # weil tlist schon von n abhängt braucht man da nur t als Index.
                        zwischenresult = [[]] * self.dim # expectation values für ein t zu jedem Level
                        dm_t[t+n*anzahl_zeitintervalle_pro_trotterintervall] = scipy.linalg.expm(L[n]*(tlist[t] - tlist[0])) @ np.reshape(dm_t[n*anzahl_zeitintervalle_pro_trotterintervall], (self.dim**2,1)) 
                        for m in range(self.dim): # für jedes t werden die expectation values zu jedem Level ausgerechnet
                            zwischenresult[m] = qutip.expect(qutip.Qobj(proj[m]), qutip.Qobj(np.reshape(dm_t[t+n*anzahl_zeitintervalle_pro_trotterintervall], (self.dim,self.dim))))
                        result[t+n*anzahl_zeitintervalle_pro_trotterintervall] = np.array(zwischenresult, dtype=np.complex128)      # für jedes t werden die 7 (je nach levelanzahl) expectation values hier eingesetzt.
                                                                                                                                    # D.h. result[t] gibt eine Liste mit den Expectation values aller Level
            result = np.array(result, dtype=np.complex128)
            if plot_pop == True:
                plot_population(tlist_ges[:], result[:], self.dim)
            print('Pop von Level {} beim Ende der Zeit {}ns, ist {}'.format(4, tlist[-1]/41341373, np.real(result[-1][3])))
            print('Pop von Level {} beim Ende der Zeit {}ns, ist {}'.format(3, tlist[-1]/41341373, np.real(result[-1][2])))
            print('Pop von Level {} beim Ende der Zeit {}ns, ist {}'.format(2, tlist[-1]/41341373, np.real(result[-1][1])))
            print('Pop von Level {} beim Ende der Zeit {}ns, ist {}'.format(7, tlist[-1]/41341373, np.real(result[-1][6])))
            print([np.real(np.amax(result[:,read_out_level])), delta_stark_shift])
            # return [np.real(np.amax(result[:,read_out_level])), delta_stark_shift]
            return [np.real((result[-1,read_out_level])), delta_stark_shift]

        elif Diagonalization == True: # Quelle http://web.mit.edu/jianshucaogroup/pdfdir/JCP117-3822-2002.pdf und https://aip.scitation.org/doi/pdf/10.1063/1.5115323 https://arxiv.org/pdf/2006.01837.pdf für die Bedingungen in equation (64)
            L = qutip.liouvillian(H, c_ops)
            L = np.reshape(L, (self.dim**2, self.dim**2))
            # # Vollständige Zeitentwicklung:
            # tlist = list(range(0,maxtime, int(maxtime/250))) # Hier wird bestimmt wie viele Werte ausgerechnet werden
            tlist = np.linspace(0, maxtime, 250)
            dm_t = [[]] * len(tlist) # Eine Liste, die für jeden Zeitpunkt die zugehörige Dichtematrix enthalten wird
            result = [[]] * len(tlist)
            # Für t=0:
            zwischenresult0 = [[]] * self.dim
            dm_t[0] = np.reshape(initial_state_dm, (self.dim**2,1))
            # dm_t[0] = np.reshape(proj[initial_state_index], (self.dim**2,1))
            # print(proj[initial_state_index])

            for m in range(self.dim): # für t=0 werden die expectation values zu jedem Level ausgerechnet
                zwischenresult0[m] = qutip.expect(qutip.Qobj(proj[m]), qutip.Qobj(np.reshape(dm_t[0], (self.dim,self.dim))))
            result[0] = np.array(zwischenresult0)
            # Für t>0:
            for t in range(1,len(tlist)): # Loopt über die Anzahl an Zeitpunkten.     Jeder Zeitpunkt wird mit t geindext
                zwischenresult = [[]] * self.dim # expectation values für ein t zu jedem Level
                dm_t[t] = scipy.linalg.expm(L*tlist[t]) @ np.reshape(dm_t[0], (self.dim**2,1)) 
                for m in range(self.dim): # für jedes t werden die expectation values zu jedem Level ausgerechnet
                    zwischenresult[m] = qutip.expect(qutip.Qobj(proj[m]), qutip.Qobj(np.reshape(dm_t[t], (self.dim,self.dim)))) # Das dürfte nicht stimmen, weil man ja in der L-Basis ist.
                result[t] = np.array(np.reshape(zwischenresult, (self.dim,)))   # für jedes t werden die 7 (je nach levelanzahl) expectation values hier eingesetzt.
                                                                                # D.h. result[t] gibt eine Liste mit den Expectation values aller Level
            result = np.array(result)
            if plot_pop == True:
                plot_population_diagonalization(tlist, result, self.dim)
            maximum_index = np.argmax(result[:,0]) # gibt Zeitpunkt, bei dem das Maximum von Level 1 erreicht wird.
            print('Pop von Level {} beim Maximum von Level 1, Zeit {}ns, ist {}'.format(4, tlist[maximum_index]/41341373, np.real(result[maximum_index][3])))  # Pop von Level 4, um zu schauen wie weit der Decay dort fortgeschritten ist
            print('Pop von Level {} beim Maximum von Level 1, Zeit {}ns, ist {}'.format(3, tlist[maximum_index]/41341373, np.real(result[maximum_index][2])))
            print('Pop von Level {} beim Maximum von Level 1, Zeit {}ns, ist {}'.format(2, tlist[maximum_index]/41341373, np.real(result[maximum_index][1])))
            print('Pop von Level {} beim Maximum von Level 1, Zeit {}ns, ist {}'.format(7, tlist[maximum_index]/41341373, np.real(result[maximum_index][6])))
            return [np.real(np.amax(result[:,read_out_level])), delta_stark_shift]

        else: # Integrations-solver von QuTip
            tlist = list(range(0,maxtime,int(41300000/100)))
            print('Länge von tlist: {}'.format(len(tlist)))
            # opts=Options(nsteps=1000)
            # print(options)
            start = time.time()
            result = qutip.mesolve(H, initial_state_dm, tlist, c_ops, proj)
            end = time.time()
            print(f'{end-start:5.3f}s')
            plot_population(result, self.dim)
            print(result.expect[1][-1])


#%% Test
import arc
from arc import *

eV_to_au = 3.67493e-02
Hz_to_au = 1.51983e-16 # f in Hz to E (bzw. omega) in au
c_au = 1.37036e02
pers_to_au = 2.41888e-17

Sr = Strontium88()
level1 = Level([(Sr.getEnergy(n=5,l=0,j=0,s=0) - Sr.getEnergy(n=5,l=0,j=0,s=0))* eV_to_au]) # 5s^2
level2 = Level([(Sr.getEnergy(n=5,l=1,j=0,s=1) - Sr.getEnergy(n=5,l=0,j=0,s=0))* eV_to_au]) # 5s5p 3P0
level3 = Level([(Sr.getEnergy(n=5,l=1,j=1,s=1) - Sr.getEnergy(n=5,l=0,j=0,s=0))* eV_to_au]) # 5s5p 3P1
level4 = Level([(Sr.getEnergy(n=5,l=1,j=2,s=1) - Sr.getEnergy(n=5,l=0,j=0,s=0))* eV_to_au]) # 5s5p 3P2
level5 = Level([(2.49824373)* eV_to_au]) # 5s4d 1D2 aus NIST
level6 = Level([(Sr.getEnergy(n=5,l=1,j=1,s=0) - Sr.getEnergy(n=5,l=0,j=0,s=0))* eV_to_au]) # 5s5p 1P1
level7 = Level([(Sr.getEnergy(n=6,l=0,j=1,s=1) - Sr.getEnergy(n=5,l=0,j=0,s=0))* eV_to_au]) # 5s6s 3S1
decay = Decay([2.01e+08*pers_to_au, 8.9e+06*pers_to_au, 2.7e+07*pers_to_au, 4.69e+04*pers_to_au, 4.2e+07*pers_to_au,\
                3.9e+03*pers_to_au, 9.4e+02*pers_to_au, 1.9e+03*pers_to_au, 9.425e-4*pers_to_au],\
                [[level6,level1], [level7,level2],[level7,level3],[level3,level1],[level7,level4],\
                [level6,level5], [level5,level4], [level5,level3], [level4,level1]])
# Parameter
Delta2 = 400e9*Hz_to_au
delta = 0
Delta1 = 0*Hz_to_au
Omega2 = 110e6*Hz_to_au
Omega3 = 110e6*Hz_to_au

# print(rabifreq(dipolelementsq(decay,3), I_0(1000e-3*watt_to_au, 350000*nm_to_au))/Hz_to_au)  # 1-3
# print(rabifreq(dipolelementsq(decay,2), I_0(100000e-6*watt_to_au, 350000*nm_to_au))/Hz_to_au)  # 3-7
# print(rabifreq(dipolelementsq(decay,1), I_0(325e-3*watt_to_au, 350000*nm_to_au))/Hz_to_au)  # 2-7

# System
laser2 = Laser(Omega2, 0.06707498275023156 + Delta2, [level2,level7])
laser3 = Laser(Omega3, 0.06622373681783161 + (Delta2-Delta1), [level3,level7])
system = System([level1, level2, level3, level4, level5, level6, level7], [laser2,laser3], decay)
system.draw()
system.fidelity([0,1,0,0,0,0,0], 0, int(1000e3*4e7), Diagonalization=True, delta_stark_shift=delta, plot_pop=True)