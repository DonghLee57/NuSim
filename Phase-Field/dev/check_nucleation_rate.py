# Silicon nucleation rate based on classical nucleation theory
# Interface energy & Nucleation rate
# Exp: Acta Materialia 54, 3227 (2006), Acta mater. 49, 439 (2001), J. Appl. Phys. 79, 15 (1996), Phys.Rev. Lett. 60, 2519 (1988),
#      Phys. Stat. Sol. (a) 4S, 313 (1978), J. Appl. Phys. 62, 1675 (1987), J. Appl. Phys. 84, 5383 (1998)
# Sim: Journal of Materials Research, 31, 3649 (2016)

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
kb_eV = constants.k/constants.e # eV/K

T_melt = 1685     # K
Vm = 1.205883199E-5 # m3/mol
# Heat of fusion ~ 50.21 kJ/mol
delta_Hf = 50500/constants.e/Vm  # eV/m3
# gamma ~ (0.35~0.45) J/m^2
gamma0 = 0.38
gamma  = gamma0/constants.e # eV/m2


def delta_G(T):
    return delta_Hf*(1-T/T_melt)

def delta_Gstar(T, gamma):
    return (16*np.pi*gamma**3)/(3*delta_G(T)**2)

def I_nuc(T,gamma):
    I0 = 1E+39 # m-3s-1
    return I0*np.exp(-(delta_Gstar(T,gamma))/(kb_eV*T))

def I_nuc2(T):
    I0 = 1.7E+44 # m-3s-1
    return I0*np.exp(-(5.3)/(kb_eV*T))

def I_nuc3(T,gamma):
    I0 = 1E+53 # m-3s-1
    return I0*np.exp(-(5.7+delta_Gstar(T,gamma))/(kb_eV*T))


fig, ax = plt.subplots(1, 4, figsize=(20,5))
T = np.linspace(700, T_melt-200)
# delta_G
ax[0].plot(T, delta_G(T))
ax[0].set_ylabel(r'$\Delta$G (eV/m3)', fontsize=12)
ax[0].set_xlabel('Temperature (K)', fontsize=12)

# delta_Gstar
ax[1].plot(T, delta_Gstar(T, gamma))
ax[1].set_ylabel(r'$\Delta$G$^*$ (eV)', fontsize=12)
ax[1].set_xlabel('Temperature (K)', fontsize=12)

# critical radius
ax[2].plot(T, 2*gamma/delta_G(T)*1E9)
ax[2].set_ylabel(r'Critical size (nm)', fontsize=12)
ax[2].set_xlabel('Temperature (K)', fontsize=12)

# Nucleation rate
# Degree celsius
if False:
    ax[3].semilogy(T-273, I_nuc(T, gamma),c='C0',ls='--')
    ax[3].semilogy(T-273, I_nuc(T, gamma/gamma0*0.3),c='C1',ls='--')
    ax[3].semilogy(T-273, I_nuc(T, gamma/gamma0*0.5),c='C2',ls='--')
    ax[3].semilogy(T-273, I_nuc2(T),c='C4',ls='--')
    ax[3].semilogy(T-273, I_nuc3(T, gamma/gamma0*0.35),c='k')
    ax[3].set_xlabel('Temperature (degree)', fontsize=12)
# Kelvin
else:
    ax[3].semilogy(T, I_nuc(T, gamma),c='C0',ls='--')
    ax[3].semilogy(T, I_nuc(T, gamma/gamma0*0.3),c='C1',ls='--')
    ax[3].semilogy(T, I_nuc(T, gamma/gamma0*0.5),c='C2',ls='--')
    ax[3].semilogy(T, I_nuc2(T),c='C4', ls='--')
    ax[3].semilogy(T, I_nuc3(T, gamma/gamma0*0.33),c='k')
    ax[3].set_xlabel('Temperature (K)', fontsize=12)
ax[3].set_ylim(1E5,1E40)
ax[3].set_ylabel(r'Nucleation rate (events/m$^3$ s)', fontsize=12)

plt.tight_layout()
plt.show()
