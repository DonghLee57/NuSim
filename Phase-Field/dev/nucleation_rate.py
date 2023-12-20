import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
kb = constants.k/constants.e # eV/K
print(kb)

T_melt = 889     # K
delta_Hf = 2.616 # eV/nm3
gamma = 0.5 # eV/nm2
def delta_G(T):
    return delta_Hf*(1-T/T_melt)

def delta_Gstar(T):
    return 16*np.pi*gamma**3/3*delta_G(T)**2

def I_nuc(T, Ea=3.3):
    I0 = 8*(1E+40) # nm-3s-1
    #print(np.exp(-(Ea + delta_Gstar(T))/(kb*T)))
    # Paper - IEEE ELECTRON DEVICE LETTERS 34 (2013)
    #return I0*np.exp(-(Ea + delta_Gstar(T))/(kb*T))
    return I0*np.exp(-(Ea + delta_Gstar(T))/(kb*T))
    #return np.exp(-2.3/kb/T)*(1-np.exp(-(delta_G(T))/kb/T))
   # return I0*np.exp(-(Ea+delta_Gstar(T))/(kb*T))*(1-np.exp(-(delta_G(T))/(kb*T)))


fig, ax = plt.subplots(1,3)
#T = np.linspace(400,T_melt+100)
T = np.linspace(300,T_melt+400)
x = 1/kb/T
ax[0].plot(T, delta_G(T))
ax[0].set_ylabel(r'$\Delta$G (eV)')
ax[1].plot(T, delta_Gstar(T))
ax[1].set_ylabel(r'$\Delta$G$^*$ (eV)')
#ax[2].semilogy(T, I_nuc(T))
#ax[2].axvline(T_melt)
ax[2].semilogy(x, I_nuc(T),c='k')
ax[2].semilogy(x, I_nuc(T, 2.3),c='b')
ax[2].set_xlim(10,35)
ax[2].set_ylim(bottom=1E-15)
ax[2].axvline(1/kb/T_melt,ls=':')
plt.tight_layout()
plt.show()
