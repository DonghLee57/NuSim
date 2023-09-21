import numpy as np
import matplotlib.pyplot as plt

MAXIT =  5001
NSAVE =  1000
Nx, Ny = 50, 50
Mx, My = Nx+2, Ny+2
dx, dy = 1, 1
dt = 0.1
L = 2.5
kappa = 0.1
W = 0.25

def g(eta,W=W):
    return W*(eta)**2*(1-eta)**2

def deriv_g(eta,W=W):
    return 2*W*eta*(1-eta)*(1-2*eta)

def LAP_2D(Grid, x, y, dx, dy, Nx, Ny):
    res = (Grid[x-1][y] + Grid[x+1][y] - 2*Grid[x][y])/dx**2\
        + (Grid[x][y-1] + Grid[x][y+1] - 2*Grid[x][y])/dy**2
    return res

if True:
    # Initializaition
    ETAS = np.zeros((Mx, My))
    R = Mx//4
    for x in range(1,Nx+1):
        for y in range(1,Ny+1):
            if (Mx//2 - x)**2 + (My//2 - y)**2 < R**2:
                ETAS[x][y] = 1.0
                
    # Evolve
    STEP = np.arange(0,MAXIT,NSAVE)
    NIMG = len(STEP)
    nid = 0
    fig, ax = plt.subplots(1, NIMG, figsize=(5*NIMG,5))
    for k in range(MAXIT):
        if k % NSAVE == 0:
            print(f"{k:d}/{MAXIT-1}")
            ax[nid].imshow(ETAS, vmin=0, vmax=1)
            ax[nid].set_title(f'Time = {dt*nid*NSAVE:.1f} unit time')
            nid += 1
        if k == MAXIT-1: break
        ETAS_new = np.zeros((Mx, My))
        Laplacian_ETAS = 0
        for x in range(1,Nx+1):
            for y in range(1,Ny+1):
                Laplacian_ETAS = LAP_2D(ETAS, x, y, dx, dy, Nx, Ny)
                ETAS_new[x][y] = ETAS[x][y] - dt*L*(deriv_g(ETAS[x][y]) - kappa*Laplacian_ETAS)
        ETAS = ETAS_new.copy()
        
    plt.tight_layout()
    plt.show()
