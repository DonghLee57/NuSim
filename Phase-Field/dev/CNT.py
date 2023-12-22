import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
kb_eV = constants.k/constants.e # eV/K

T_melt= 1685 # K
gamma = 0.45/constants.e/1E+18 # eV/nm^2
# J/m^2 -> eV/m^2 -> eV/nm^2
#2.81 # eV/nm^2
dHc   = 0.12 # eV
dHc   = 11.9/constants.e/constants.N_A # eV/nm^3
#dHc   = 11.9/constants.e/constants.N_A/2.555*1E+3 # eV/nm^3
# J/mol -> eV/mol -> eV/atom -> eV/angs3 -> eV/nm^3
print(dHc)

def dG(T):
    return dHc*(1-T/T_melt)

def dGstar(T):
    return 16*np.pi*gamma**3/(3*dG(T)**2)

def r_crt(T):
    return 2*gamma/dG(T)

fig, ax = plt.subplots(1,2, figsize=(10,5))
T = np.linspace(500,1500)
ax[0].plot(T, dGstar(T))
#print(f'RTA: Linear size ~ {(growth(T)/rate_nuc(T))**(0.25)*1E+9:.2f} nm')
ax[0].set_ylabel('Nucleation barrier (eV)')
#ax[0].set_ylim([0, 1])
ax[0].set_xlabel('Temperature (K)')
#ax[0].set_xlim(left=0)
#ax[0].legend()

ax[1].plot(T, r_crt(T))
ax[1].set_ylabel('Critical radius (nm)')
ax[1].set_xlabel('Temperature (K)')
#ax[1].set_ylim([0, 1])
#ax[1].set_xlabel('time (hour)')
#ax[1].set_xlim(left=0)
#ax[1].legend()
plt.tight_layout()
plt.show()
