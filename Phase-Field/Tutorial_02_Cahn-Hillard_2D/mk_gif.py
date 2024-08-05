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

fig,ax = plt.subplots(1,2)
ax[0].plot(free_energy)
ax[1].plot(conc)
ax[1].set_ylim([0.48, 0.52])

plt.show()
