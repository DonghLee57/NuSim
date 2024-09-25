from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm
from fipy.tools import numerix
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Mesh grid setting
nx = 50
dx = 1.
L  = nx*dx
mesh = Grid1D(dx=dx, nx=nx)

# Initialization
phi = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = 0.)
x = mesh.cellCenters[0]
phi.setValue(1, where=(x < 0.5*L))

D = 1.

# Governing equation
eq = TransientTerm() == DiffusionTerm(coeff=D)

timeStepDuration = 0.9 * dx**2 / (2 * D)

steps = 1000
vstep = 100
data = np.zeros((steps//vstep, nx))
total = np.zeros(steps//vstep)
for step in range(steps):
    eq.solve(var=phi,
             dt=timeStepDuration)
    if step % vstep == 9: 
        X = numerix.array(x)
        Y = numerix.array(phi)
        data[step//vstep] = Y
        total[step//vstep] = integrate.simps(Y, x=X)
        
fig, ax = plt.subplots(1, 2, figsize=(10,5), constrained_layout=True)
for i in range(steps//vstep):
    ax[0].plot(X, data[i])
ax[1].scatter(np.arange(steps//vstep), total)
ax[1].set_ylim([20,30])
plt.show()
