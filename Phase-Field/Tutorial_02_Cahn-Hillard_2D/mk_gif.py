import pyvista as pv
import numpy as np

load_file = "simulation_results.pvd"
save_file = "res.mp4"

def grad(arr):
    return (np.roll(arr, -1, axis=0) - np.roll(arr, 1, axis=0))/2

pvd = pv.get_reader(load_file)
sargs = dict(height=0.40, vertical=True, position_x= 0.8, position_y=0.4, n_labels=5, title='psi', fmt='%.2f')
mesh_kwargs = {'scalars':'psi',
              'scalar_bar_args':sargs,
              'clim':[0,1],
              'cmap':'bwr'}

p = pv.Plotter()
fformat = save_file.split('.')[-1]
if fformat == 'gif':
    # pip install imageio
    p.open_gif(save_file)
elif fformat == 'mp4':
    # pip install imageio[ffmpeg]
    # pip install imageio[pyav]
    p.open_movie(save_file, framerate=16)

nframe = len(pvd.time_values)
free_energy = []
conc = []
for frame in range(nframe):
    data = pvd.set_active_time_point(frame)
    grid = pvd.read()[0]
    arr = grid.point_data.get_array(mesh_kwargs['scalars'])

    grad_x = grad(arr)
    grad_y = grad(arr.T).T
    grad_mag = np.sqrt(grad_x**2 + grad_y**2)
    free_energy.append(np.sum( arr**2 * (1-arr)**2 + 0.5/2 * grad_mag))
    conc.append(np.sum(arr)/(64**2))

    p.add_mesh(pvd.read()[0], show_scalar_bar=True, **mesh_kwargs)
    p.camera_position = 'xy'
    p.write_frame()

p.close()
    
import matplotlib.pyplot as plt
fs = 12
fig,ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
x = np.linspace(0, 1, len(free_energy))
ax[0].plot(x, free_energy)
ticks = np.arange(0.0, 1.1, 0.25)
ax[0].set_xticks(ticks)
ax[0].set_xticklabels([f"{x:.1f}" for x in ticks], fontsize=fs)
ax[0].set_xlim([0, 1])
ax[0].set_xlabel("Time", fontsize=fs)
ticks = np.arange(0.0, 401, 100)
ax[0].set_yticks(ticks)
ax[0].set_yticklabels([f"{x:.0f}" for x in ticks], fontsize=fs)
ax[0].set_ylim([0, 400])
ax[0].set_ylabel("Free energy", fontsize=fs)


ax[1].plot(x, conc)
ticks = np.arange(0.0, 1.1, 0.25)
ax[1].set_xticks(ticks)
ax[1].set_xticklabels([f"{x:.1f}" for x in ticks], fontsize=fs)
ax[1].set_xlim([0, 1])
ax[1].set_xlabel("Time", fontsize=fs)
ticks = np.arange(0.48, 0.521, 0.01)
ax[1].set_yticks(ticks)
ax[1].set_yticklabels([f"{x:.2f}" for x in ticks], fontsize=fs)
ax[1].set_ylim([0.48, 0.52])
ax[1].set_ylabel(r"Average psi per grid point", fontsize=fs)

plt.show()
