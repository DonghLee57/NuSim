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

def I_nuc(T, Ea=3.3):
    I0 = 8*(1E+40) # nm-3s-1
    return I0*np.exp(-(Ea + delta_Gstar(T))/(kb_eV*T))

fig, ax = plt.subplots(1,3, figsize=(15,5))
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
plt.tight_layout()
plt.show()
