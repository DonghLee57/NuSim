import pyvista as pv

pvd = pv.get_reader("simulation_results.pvd")
sargs = dict(height=0.40, vertical=True, position_x= 0.8, position_y=0.4, n_labels=5, title='psi', fmt='%.2f')
mesh_kwargs = {'scalars':'psi',
              'scalar_bar_args':sargs,
              'clim':[0,1],
              'cmap':'bwr'}

p = pv.Plotter()
fname = 'res.mp4'
fformat = fname.split('.')[-1]
if fformat == 'gif':
    # pip install imageio
    p.open_gif("res.gif")
elif fformat == 'mp4':
    # pip install imageio[ffmpeg]
    # pip install imageio[pyav]
    p.open_movie("res.mp4", framerate=16)

nframe = len(pvd.time_values)
for frame in range(nframe):
    data = pvd.set_active_time_point(frame)
    p.add_mesh(pvd.read()[0], show_scalar_bar=True, **mesh_kwargs)
    p.camera_position = 'xy'
    p.write_frame()

p.close()
    
