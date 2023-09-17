import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

Ngr = int(sys.argv[1])
Nx, Ny = 100, 100
#

IMGS =  np.arange(0, 301, 100)
fig, ax= plt.subplots(1,len(IMGS), figsize=(5*len(IMGS), 5))
for k in range(len(IMGS)):
    SHOW = np.zeros((Nx, Ny))
    SHOW_p = np.zeros((Nx, Ny))
    for n in range(Ngr):
        FILE = f'PhaseField_{IMGS[k]:d}_{n:d}.dat'
        with open(FILE, 'r') as o: tmp = o.readlines()
        Grid = []
        for i in range(len(tmp)):
            Grid.append(list(map(float,tmp[i].split())))
        Grid = np.array(Grid)
        Grid[Grid[:,:] < 0.8] = 0
        SHOW += Grid
         
    ax[k].imshow(SHOW, vmin=0, vmax=1)
    ax[k].set_ylabel('cpp')

plt.show()
