# Cahn-Hilliard equation in phae-field modeling with finite difference algorithm

import numpy as np
import matplotlib.pyplot as plt

#
MAXIT = 5001
Nx, Ny = 64, 64       
dx, dy = 1.0, 1.0     
dt = 0.01

#
Conc0 = 0.4
noise = 0.02
Mob   = 1.0
grad_coeff = 0.5
A = 1.0
plots = ['static']

# Initializaition
Conc = np.zeros((MAXIT,Nx,Ny))
for x in range(Nx):
  for y in range(Ny):
    Conc[0][x][y] = Conc0 + noise*(0.5 - np.random.random())

# Evolve
for k in range(MAXIT-1):
  if k % 100 ==0: print(f"{k+1:d}/{MAXIT-1}")
  Laplacian_C = np.zeros((Nx,Ny))
  Laplacian_dFdc = np.zeros((Nx,Ny))
  dfdc = np.zeros((Nx,Ny))
  dFdc = np.zeros((Nx,Ny))
  for x in range(Nx):
    for y in range(Ny):
      idx = (x+1) % Nx
      idy = (y+1) % Ny
      Laplacian_C[x][y] = (Conc[k][idx-2][y] + Conc[k][idx][y] + Conc[k][x][idy-2] + Conc[k][x][idy] - 4*Conc[k][x][y])/(dx*dy)
      dfdc[x][y] = 2*A*Conc[k][x][y]*(1-Conc[k][x][y])*(1-2*Conc[k][x][y])
      dFdc[x][y] = dfdc[x][y] - grad_coeff*Laplacian_C[x][y]
  for x in range(Nx):
    for y in range(Ny):
      idx = (x+1) % Nx
      idy = (y+1) % Ny
      Laplacian_dFdc[x][y] = (dFdc[idx-2][y] + dFdc[idx][y] + dFdc[x][idy-2] + dFdc[x][idy] - 4*dFdc[x][y])/(dx*dy)
      Conc[k+1][x][y] = Conc[k][x][y] + dt*Mob*Laplacian_dFdc[x][y]

if 'static' in plots:
  NIMG = 11
  fig, ax = plt.subplots(1,NIMG, figsize=((NIMG)*5,5))
  STEP = MAXIT//NIMG#//10
  for n in range(NIMG):
    ax[n].imshow(Conc[n*STEP])
    ax[n].set_title(f'Time = {dt*n*STEP:.1f} unit time')
  plt.show()
