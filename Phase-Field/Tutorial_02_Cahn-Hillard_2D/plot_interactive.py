import pyvista as pv

pvd = pv.get_reader("simulation_results.pvd")
sargs = dict(height=0.40, vertical=True, position_x= 0.8, position_y=0.4, n_labels=5, title='psi', fmt='%.2f')
mesh_kwargs = {'scalars':'phase_field',
              'scalar_bar_args':sargs,
              'clim':[0,1],
              'cmap':'bwr'}

def show_grid(step):
    data = pvd.set_active_time_point(int(step)-1)
    p.add_mesh(pvd.read()[0], **mesh_kwargs)
    p.camera_position = 'xy'

p = pv.Plotter()
p.add_slider_widget(show_grid, rng=(1,25), title='psi')
p.show()
