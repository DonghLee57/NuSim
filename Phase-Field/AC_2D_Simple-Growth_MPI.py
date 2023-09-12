# Allenn-Cahn equation in phae-field modeling with finite difference algorithm
#

import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#
MAXIT = 1001
NSAVE = 100
Nx, Ny = 64, 64
dx, dy = 0.5, 0.5
dt = 0.005

#
Mob   = 5.0
grad_coeff = 0.1
A, B = 1.0, 1.0

def main():
    # Initializaition
    Ngr = 1
    ETAS = np.zeros((Ngr, Nx, Ny))
    R = 5
    np.random.seed(1234)
    if rank == 0 :
        for N in range(Ngr):
            x0 = np.random.randint(low=0,high=Nx)
            y0 = np.random.randint(low=0,high=Ny)
            for x in range(Nx):
                for y in range(Ny):
                    # Periodic boundary
                    if np.fabs(x0-x) > Nx//2: x += np.sign(x0-x)*Nx
                    if np.fabs(y0-y) > Ny//2: y += np.sign(y0-y)*Ny
                    if ((x0 - x)**2 + (y0 - y)**2)**0.5 < R:
                        ETAS[N][x%Nx][y%Ny] = 1.0
    comm.Barrier()
    ETAS = comm.bcast(ETAS)

    # Evolve
    for k in range(MAXIT):
        ETAS_new = np.zeros((Ngr, Nx, Ny))
        if k % NSAVE == 0: 
            if rank == 0:
                print(f"{k:d}/{MAXIT-1}")
                with open(f"data{k:05d}.pickle",'wb') as f:
                    pickle.dump(ETAS, f)
        if k == MAXIT-1: break
        for N in range(Ngr):
            Laplacian_ETAS = np.zeros((Nx,Ny))
            df = np.zeros((Nx,Ny))

            psize = Nx // size
            if rank == size - 1:
                for x in range(rank*psize,Nx):
                    for y in range(Ny):
                        idx = (x+1) % Nx
                        idy = (y+1) % Ny
                        Laplacian_ETAS[x][y] = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny)
                        SUM_eta_sq = 0
                        for ngr in range(Ngr):
                            if N != ngr: SUM_eta_sq += ETAS[ngr][x][y]**2
                        df[x][y] = -A*ETAS[N][x][y] + B*ETAS[N][x][y]**3 + 2*ETAS[N][x][y]*SUM_eta_sq
            else:
                for x in range(rank*psize,(rank+1)*psize):
                    for y in range(Ny):
                        idx = (x+1) % Nx
                        idy = (y+1) % Ny
                        Laplacian_ETAS[x][y] = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny)
                        SUM_eta_sq = 0
                        for ngr in range(Ngr):
                            if N != ngr: SUM_eta_sq += ETAS[ngr][x][y]**2
                        df[x][y] = -A*ETAS[N][x][y] + B*ETAS[N][x][y]**3 + 2*ETAS[N][x][y]*SUM_eta_sq

            if rank == size - 1:
                for x in range(rank*psize,Nx):
                    for y in range(Ny):
                        ETAS_new[N][x][y] = ETAS[N][x][y] - dt*Mob*(df[x][y] - grad_coeff*Laplacian_ETAS[x][y])
            else:
                for x in range(rank*psize,(rank+1)*psize):
                    for y in range(Ny):
                        ETAS_new[N][x][y] = ETAS[N][x][y] - dt*Mob*(df[x][y] - grad_coeff*Laplacian_ETAS[x][y])
      
        ETAS_new = comm.reduce(ETAS_new)
        comm.Barrier()
        ETAS_new = comm.bcast(ETAS_new)
        ETAS = ETAS_new.copy()
    return 0

def lap_2D(Grid, x, y, dx, dy, Nx, Ny):
    idx = (x+1) % Nx
    idy = (y+1) % Ny
    res = (Grid[idx-2][y] + Grid[idx][y] - 2*Grid[x][y])/dx**2\
        + (Grid[x][idy-2] + Grid[x][idy] - 2*Grid[x][y])/(dy**2)
    return res

def plots(PATH='./', NIMG=11):
    fig, ax = plt.subplots(1,NIMG, figsize=((NIMG)*5,5))
    STEP = MAXIT//NIMG
    for n in range(NIMG):
        SHOW = np.zeros((Nx,Ny))
        with open(f'{PATH}/data{n*NSAVE:05d}.pickle','rb') as f:
            Grid = pickle.load(f)
        for N in range(len(Grid)):
            SHOW += Grid[N][:]*(N+1)
        ax[n].imshow(SHOW)
        ax[n].set_title(f'Time = {dt*n*STEP:.1f} unit time')
    plt.show()
    return 0

if __name__ == "__main__":
    main()
    plots()
