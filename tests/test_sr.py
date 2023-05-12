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