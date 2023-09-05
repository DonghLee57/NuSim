import numpy as np
import matplotlib.pyplot as plt

# rho: density
# cp: heat capacity
# lambda: thermal conductivity
# mu = lambda/rho/cp

# 
MAXIT = 10000
Nx = 1000
plots = ['static', 'ani'] 
#

TEMP = np.empty((MAXIT, Nx))

# mu*dt/dx**2 < 1/2
mu = 2.0  # mm2/s
dt = 0.2  # s
dx = 1    # mm
alpha = mu*dt/dx**2
if stability <= 0.5:
  print(f"Simulation condition is fine. Stability is {alpha:.2f} under 0.5.")
else:
  print(f"Simulation condition is bad.  Stability is {alpha:.2f}  over 0.5.")

# Initializaition
TEMP[0][Nx*4//10:Nx*6//10] = 1000

# Calculate Heat Eq.
for k in range(MAXIT-1):
  for x in range(Nx):
    # "idx" is defined for Periodic Boundary Condition
    idx = (x+1) % Nx
    TEMP[k+1][idx-1] = TEMP[k][idx-1] + alpha * (TEMP[k][idx] - 2*TEMP[k][idx-1] + TEMP[k][idx-2])

if 'static' in plots:
  NIMG = 10
  fig, ax = plt.subplots(1,NIMG, figsize=(NIMG*5,5))
  STEP = MAXIT//NIMG
  for n in range(NIMG):
    plt_T = ax[n].plot(np.arange(0,dx*Nx, dx)/1000, TEMP[n*STEP])
    ax[n].set_title(f'Time = {dt*n:.1f} sec')
    ax[n].set_xlabel('Position (m)')
    ax[n].set_xlabel('Temperature (K)')
    ax[n].set_ylim([-10,1010])
  plt.show()

if 'ani' in plots:
  from matplotlib import animation
  NIMG = 10
  STEP = MAXIT//NIMG
  fig, ax = plt.subplots(figsize=(5,5))
  Writer = animation.writers['ffmpeg']
  metadata = dict(title='Movie', artist='Matplotlib')
  writer = Writer(fps=1, metadata=metadata)
  with writer.saving(fig, "Heat_1D.mp4", 100):
    for n in range(NIMG):
      plt_T = ax.plot(np.arange(0,dx*Nx, dx)/1000, TEMP[n*STEP])
      ax.set_title(f'Time = {dt*n:.1f} sec')
      ax.set_xlabel('Position (m)')
      ax.set_xlabel('Temperature (K)')
      ax.set_ylim([-10,1010])
      writer.grab_frame()
      plt.cla()
