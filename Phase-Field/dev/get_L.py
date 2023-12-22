# Journal of Crystal Growth 585, 126590 (2022), ACS Omega 8, 3329 (2023)

# gamma: interfacial energy
# delta: interfacial width
# b: ~2.2

def D(T):
    D0 = 1E-4 # cm2/s
    Ea = 2.0  # eV
    return D0*np.exp(-Ea/(kb_eV*T))

def A(T):
    # Relation between A and D
    # A = -2.0E-22*D**(-0.78)
    return -2.0E-22*D(T)**(-0.78) 

def L(T):
    # Unit of L: distance^3/energy/time
    #A = 1 # m^4/(J*s)
    B = 1 # unitless (0.6~1.4)
    T_melt = 1000
    delta_Tr = (1-T/T_melt)
    spacing = 2.4E-10 # m
    return A(T)*(D(T)*delta_Tr/(2*np.pi*spacing**2))**B

def K():
    return 12*delta*gamma/b

def W():
    return 6*gamma*b/delta

