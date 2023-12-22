# Journal of Crystal Growth 585, 126590 (2022), ACS Omega 8, 3329 (2023)

# gamma: interfacial energy
# delta: interfacial width
# b: ~2.2

def L(T):
    A = 1 # m^4/(J*s)
    B = 1 # unitless
    T_melt = 1000
    delta_Tr = (1-T/T_melt)
    spacing = 2.4 # angstrom
    return A*(D(T)*delta_Tr/(2*np.pi*spacing**2))**B

def K():
    return 12*delta*gamma/b

def W():
    return 6*gamma*b/delta

def D(T):
    D0 = 1E-4 # cm2/s
    Ea = 2.0  # eV
    return D0*np.exp(-Ea/(kb_eV*T))
