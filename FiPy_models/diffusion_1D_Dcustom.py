from fipy import CellVariable, FaceVariable, Grid1D, TransientTerm, DiffusionTerm
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
cx = mesh.cellCenters[0]
phi.setValue(1, where=(cx < 0.5*L))

D = FaceVariable(mesh = mesh,
                 value = 1.0)
fx = mesh.faceCenters[0]
D.setValue(0.1, where=(0.5*L < fx))

# Governing equation
eq = TransientTerm() == DiffusionTerm(coeff=D)

timeStepDuration = 0.9 * dx**2 / (2 * max(numerix.array(D)))

steps = 1000
vstep = 100
data = np.zeros((steps//vstep, nx))
total = np.zeros(steps//vstep)
for step in range(steps):
    eq.solve(var=phi,
             dt=timeStepDuration)
    if step % vstep == 9: 
        X = numerix.array(cx)
        Y = numerix.array(phi)
        data[step//vstep] = Y
        total[step//vstep] = integrate.simps(Y, x=X)

fig, ax = plt.subplots(1, 2, figsize=(10,5), constrained_layout=True)
for i in range(steps//vstep):
    ax[0].plot(X, data[i])
ax[1].scatter(np.arange(steps//vstep), total)
ax[1].set_ylim([20,30])
plt.show()
