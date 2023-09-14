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
Nx, Ny = 100, 100
dx, dy = 1, 1

# Boundary condition
# 'p' : periodic
# 'f' : fixed
x0, xN = 'f', 'f'
y0, yN = 'f', 'f'

dt = 0.002

# Mob is independent on Nx, Ny
Mob   = 10.0
grad_coeff = 0.1
A, B = 1.0, 1.0

def main():
    # Initializaition
    Ngr = 10
    ETAS = np.zeros((Ngr, Nx, Ny))
    R = 25
    #np.random.seed(1234)
    if rank == 0 :
        for N in range(Ngr):
            xc = np.random.randint(low=0,high=Nx)
            yc = np.random.randint(low=0,high=Ny)
            for x in range(Nx):
                for y in range(Ny):
                    if x0[0] == 'p':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) < 0: x -= Nx
                    elif x0[0] == 'f':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) < 0: x = xc
                    if xN[0] == 'p':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) > 0: x += Nx
                    elif xN[0] == 'f':
                        if np.fabs(xc-x) > Nx//2 and np.sign(xc-x) > 0: x = xc
                    if y0[0] == 'p':
                        if np.fabs(yc-y) > Ny//2 and np.sign(yc-y) < 0: y -= Ny
                    elif y0[0] == 'f':
                        if np.fabs(yc-y) > Nx//2 and np.sign(yc-y) < 0: y = yc
                    if yN[0] == 'p':
                        if np.fabs(yc-y) > Ny//2 and np.sign(yc-y) > 0: y += Ny
                    elif y0[0] == 'f':
                        if np.fabs(yc-y) > Nx//2 and np.sign(yc-y) > 0: y = yc
                    if ((xc - x)**2 + (yc - y)**2)**0.5 < R:
                        ETAS[N][x%Nx][y%Ny] = 1.0
    ETAS = comm.reduce(ETAS)
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
                        Laplacian_ETAS[x][y] = lap_2D(ETAS[N], x, y, dx, dy, Nx, Ny)
                        SUM_eta_sq = 0
                        for ngr in range(Ngr):
                            if N != ngr: SUM_eta_sq += ETAS[ngr][x][y]**2
                        df[x][y] = -A*ETAS[N][x][y] + B*ETAS[N][x][y]**3 + 2*ETAS[N][x][y]*SUM_eta_sq
            else:
                for x in range(rank*psize,(rank+1)*psize):
                    for y in range(Ny):
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
    if x0[0] == 'p':   idx0 = x - 1
    elif x0[0] == 'f':
        if x == 0: idx0 = x
        else:      idx0 = x - 1
    if xN[0] == 'p':   idxN = (x+1)%Nx
    elif xN[0] == 'f':
        if x == Nx-1: idxN = x
        else:         idxN = x + 1
    if y0[0] == 'p':   idy0 = y - 1
    elif y0[0] == 'f':
        if y == 0: idy0 = y
        else:      idy0 = y - 1
    if yN[0] == 'p':   idyN = (y+1)%Ny
    elif yN[0] == 'f':
        if y == Ny-1: idyN = y
        else:         idyN = y + 1
    res = (Grid[idx0][y] + Grid[idxN][y] - 2*Grid[x][y])/dx**2\
        + (Grid[x][idy0] + Grid[x][idyN] - 2*Grid[x][y])/dy**2
    return res

def plot_res(PATH='./', nrow=1):
    STEP = np.arange(0,MAXIT,NSAVE)
    NIMG = len(STEP)
    if NIMG == 1:
        fig, ax = plt.subplots(1, NIMG, figsize=(5,5))
        SHOW = np.zeros((Nx,Ny))
        with open(f'{PATH}/data{0:05d}.pickle','rb') as f:
            Grid = pickle.load(f)
        for N in range(len(Grid)):
            SHOW += Grid[N][:]*(N+1)
        ax.imshow(SHOW, vmin=0, vmax=len(Grid))
        ax.set_title(f'Time = {0:.1f} unit time')
        ax.set_xlabel(f'{dx*Nx:.2f} unit distance')
    else:
        if nrow > 1:
            if NIMG % nrow != 0:
                fig, ax = plt.subplots(nrow+1, NIMG//nrow, figsize=(NIMG//nrow*5,5*nrow))
            if NIMG % nrow == 0:
                fig, ax = plt.subplots(nrow, NIMG//nrow, figsize=(NIMG//nrow*5,5*nrow))
            for n in range(NIMG):
                SHOW = np.zeros((Nx,Ny))
                with open(f'{PATH}/data{n*NSAVE:05d}.pickle','rb') as f:
                    Grid = pickle.load(f)
                for N in range(len(Grid)):
                    SHOW += Grid[N][:]*(N+1)
                ax[n//(NIMG//nrow)][n%(NIMG//nrow)].imshow(SHOW, vmin=0, vmax=len(Grid))
                ax[n//(NIMG//nrow)][n%(NIMG//nrow)].set_title(f'Time = {dt*n*NSAVE:.1f} unit time')
                ax[n//(NIMG//nrow)][n%(NIMG//nrow)].set_xlabel(f'{dx*Nx:.2f} unit distance')
        elif nrow == 1:
            fig, ax = plt.subplots(1, NIMG, figsize=(NIMG*5,5))
            for n in range(NIMG):
                SHOW = np.zeros((Nx,Ny))
                with open(f'{PATH}/data{n*NSAVE:05d}.pickle','rb') as f:
                    Grid = pickle.load(f)
                for N in range(len(Grid)):
                    SHOW += Grid[N][:]*(N+1)
                ax[n].imshow(SHOW, vmin=0, vmax=len(Grid))
                ax[n].set_title(f'Time = {dt*n*NSAVE:.1f} unit time')
                ax[n].set_xlabel(f'{dx*Nx:.2f} unit distance')
    plt.tight_layout()
    plt.show()
    return 0

def check_stability():
    res = (Mob*dt)*(1/dx**2+1/dy**2)
    if res < 0.25:
        print(f"Maybe stable. {res:.2f} is less than 0.25.")
        return True
    else: 
        print(f"Maybe unstable. {res:.2f} is greater than 0.25.")
        return False

if __name__ == "__main__":
    if check_stability():
        main()
        if rank==0: plot_res(nrow=3)
    else:
        sys.exit()
