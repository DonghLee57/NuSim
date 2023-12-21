# Paper - IEEE ELECTRON DEVICE LETTERS 34 (2013)

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
kb_eV = constants.k/constants.e # eV/K

T_melt = 889     # K
delta_Hf = 2.616 # eV/nm3
gamma = 0.5 # eV/nm2
def delta_G(T):
    return delta_Hf*(1-T/T_melt)

def delta_Gstar(T):
    return (16*np.pi*gamma**3)/(3*delta_G(T)**2)

def I_nuc(T):
    I0 = 8E+45 # nm-3s-1
    Ea = 3.3
    return I0*np.exp(-(Ea + delta_Gstar(T))/(kb_eV*T))

def Vg(T):
    V0 = 5E+24
    Ea = 2.3
    return V0*np.exp(-Ea/(kb_eV*T))

def Xc(T,time):
    return 1-np.exp(-np.pi*Vg(T)**3*I_nuc(T)*time**4/3)

fig, ax = plt.subplots(1,4, figsize=(20,5))
T = np.linspace(300,850)
x = 1/kb_eV/T
ax[0].plot(T, delta_G(T))
ax[0].set_ylabel(r'$\Delta$G (eV)')
ax[1].plot(T, delta_Gstar(T))
ax[1].set_ylabel(r'$\Delta$G$^*$ (eV)')
ax[2].semilogy(x, I_nuc(T),c='k')
ax[2].set_xlim(10,35)
ax[2].set_ylim(bottom=1E-15)
ax[2].axvline(1/kb_eV/T_melt,ls=':',c='k')
time = np.linspace(0,50)
ax[3].plot(time, Xc(190+273, time))
ax[3].set_ylim([0,1])
plt.tight_layout()
plt.show()
