# Allenn-Cahn equation in phae-field modeling with finite difference algorithm
#

import numpy as np
import matplotlib.pyplot as plt

#
MAXIT = 1001
Nx, Ny = 64, 64
dx, dy = 0.5, 0.5
dt = 0.005

#
Mob   = 5.0
grad_coeff = 0.1
A, B = 1.0, 1.0
plots = ['static']

# Initializaition
Ngr = 3
ETA = np.zeros((MAXIT, Ngr, Nx, Ny))
ETAS = np.zeros((MAXIT, Ngr, Nx, Ny))
R = 10
for N in range(Ngr):
  x0 = np.random.randint(low=0,high=Nx)
  y0 = np.random.randint(low=0,high=Ny)
  for x in range(Nx):
    for y in range(Ny):
      # Periodic boundary
      if np.fabs(x0-x) > Nx//2: x += np.sign(x0-x)*Nx
      if np.fabs(y0-y) > Ny//2: y += np.sign(y0-y)*Ny
      if ((x0 - x)**2 + (y0 - y)**2)**0.5 < R:
          ETAS[0][N][x%Nx][y%Ny] = 1.0

# Evolve
for k in range(MAXIT-1):
  if k % 100 == 99: print(f"{k+1:d}/{MAXIT-1}")
  for N in range(Ngr):
    Laplacian_ETAS = np.zeros((Nx,Ny))
    df = np.zeros((Nx,Ny))
    for x in range(Nx):
      for y in range(Ny):
        idx = (x+1) % Nx
        idy = (y+1) % Ny
        Laplacian_ETAS[x][y] = (ETAS[k][N][idx-2][y] + ETAS[k][N][idx][y] + ETAS[k][N][x][idy-2] + ETAS[k][N][x][idy] - 4*ETAS[k][N][x][y])/(dx*dy)
        SUM_eta_sq = 0
        for ngr in range(Ngr):
          if N != ngr: SUM_eta_sq += ETAS[k][ngr][x][y]**2
        df[x][y] = -A*ETAS[k][N][x][y] + B*ETAS[k][N][x][y]**3 + 2*ETAS[k][N][x][y]*SUM_eta_sq

    for x in range(Nx):
      for y in range(Ny):
        ETAS[k+1][N][x][y] = ETAS[k][N][x][y] - dt*Mob*(df[x][y] - grad_coeff*Laplacian_ETAS[x][y])
        if ETAS[k+1][N][x][y] > 1: print(k+1,N, x,y, ETAS[k+1][N][x][y])

if 'static' in plots:
  NIMG = 11
  fig, ax = plt.subplots(1,NIMG, figsize=((NIMG)*5,5))
  STEP = MAXIT//NIMG#//10
  for n in range(NIMG):
    SHOW = np.zeros((Nx,Ny))
    for N in range(Ngr):
      SHOW += ETAS[n*STEP][N][:]*(N+1)
    ax[n].imshow(SHOW)
    ax[n].set_title(f'Time = {dt*n*STEP:.1f} unit time')
  plt.show()
