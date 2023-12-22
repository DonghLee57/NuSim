import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
kb_eV = constants.k/constants.e # eV/K
# Si molar volume
Vmolar = 1.205883E-5*1E+27 # nm^3/mol

T_melt= 1685 # K
# 340~450 mJ/m2
gamma = 0.45/constants.e/1E+18 # eV/nm^2
print(f"Interfacial energy: {gamma:.2f} eV/nm^2")
dHc   = 11900/constants.e/Vmolar # eV/nm^3 Heat of crystallization
#dHc   = 50550/constants.e/Vmolar # eV/nm^3 Heat of fusion
# J/mol -> eV/mol -> eV/nm^3
print(f"Latent heat of crystallization: {dHc:.2f} eV/nm^3")

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
ax[0].set_xlabel('Temperature (K)')
ax[0].set_ylim(bottom=0)

ax[1].plot(T, 2*r_crt(T))
ax[1].set_ylabel('Critical size (nm)')
ax[1].set_xlabel('Temperature (K)')
ax[1].set_ylim(bottom=0)

plt.tight_layout()
plt.show()
