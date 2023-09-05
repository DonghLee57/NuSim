# Non-conserved Allen-Cahn equation in phae-field modeling
#

import numpy as np
import matplotlib.pyplot as plt

def free_energy(Nx, Ny, Ngr, ETA, ETAS, ig):
  A, B = 1, 1
  SUM = np.zeros((Nx, Ny))
  for g in range(Ngr):
      if g != ig:
        SUM += ETAS[g][:,:]**2
  df = A*(2*B*ETA*SUM)+ETA**3-ETA
  return df


# eta: order paramter  (0, 1)
# L: mobility coefficient
# eps: the gradient energy coefficient
# A, B: Constants

#
MAXIT = 1001
Ngr = 1
Nx, Ny = 64, 64       #
dx, dy = 0.5, 0.5     #
dt = 0.005

#
L = 5.0
eps = 0.1  
plots = ['static']
#

ETA = np.zeros((MAXIT, Nx, Ny))
ETAS= np.zeros((Ngr, Nx, Ny))

# Initializaition
ETAS[0][:Nx//4,:Ny//4].fill(1)
for g in range(Ngr):
  ETA[0] += ETAS[g] 

# Calculate Allen-Cahn eqaution
for k in range(MAXIT-1):
  lap_ETA = np.zeros((Nx,Ny))
  for g in range(Ngr):
    for x in range(Nx):
      for y in range(Ny):
        idx = (x+1) % Nx
        idy = (y+1) % Ny
        lap_ETA[idx-1][idy-1] = (ETA[k][idx-2][idy-1] + ETA[k][idx][idy-1] + ETA[k][idx-1][idy-2] + ETA[k][idx-2][idy] - 4*ETA[k][idx-1][idy-1])/(dx*dy)
    dfdeta = free_energy(Nx, Ny, Ngr, ETA[k], ETAS, g)
    ETA[k+1] += ETA[k] - dt * L * (dfdeta - eps * lap_ETA)

if 'static' in plots:
  NIMG = 11
  fig, ax = plt.subplots(1,NIMG, figsize=((NIMG)*5,5))
  STEP = MAXIT//NIMG#//10
  for n in range(NIMG):
    ax[n].imshow(ETA[n*STEP])
    ax[n].set_title(f'Time = {dt*n:.1f} unit time')
  plt.show()
