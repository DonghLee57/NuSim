import numpy as np
import matplotlib.pyplot as plt

SIZE  = 50
MAXIT = 10

alpha = 2.0
dx = 1
dy = 1
dt = dx**2/4/alpha

DATA = np.zeros((MAXIT, SIZE, SIZE))
DATA[:, (SIZE-1):, :] = 100

for t in range(0, MAXIT-1, 1):
  for i in range(1, SIZE-1, dx):
    for j in range(1, SIZE-1, dy):
      DATA[t + 1, i, j] = dt*alpha*((DATA[t][i+1][j] - 2*DATA[t][i][j] + DATA[t][i-1][j])/dx**2 + (DATA[t][i][j+1] - 2*DATA[t][i][j] + DATA[t][i][j-1])/dy**2) + DATA[t][i][j]

"""
fig, ax = plt.subplots()
ax.set_title(f"Temperature at t = {k*delta_t:.3f} unit time")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=0, vmax=100)
plt.colorbar()
plt.clf()
"""
