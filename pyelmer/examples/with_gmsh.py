import os
import gmsh
from pyelmer import elmer, execute
from pyelmer.gmsh_utils import add_physical_group

# geometry modeling
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add('heat flux')

rect =  gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
gmsh.model.occ.synchronize()

ph_rect = add_physical_group(2, [rect], 'rect1')
bndry_left = add_physical_group(1, [4], 'left')
bndry_right = add_physical_group(1, [2], 'right')
# Here, the IDs of the boundaries (4 - left, 2 - right) were taken from
# the visualization. More advanced functions are available in
# pyelmer.gmsh_utils. 

# mesh
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.1)
gmsh.model.mesh.generate(2)
gmsh.write('first_order.msh2')

# visualize
gmsh.fltk.run()

# export
if not os.path.exists('./simdata'):
    os.mkdir('./simdata')
gmsh.write('./simdata/rectangle.msh2')

# convert
execute.run_elmer_grid('./simdata', 'rectangle.msh2')

# elmer setup
sim = elmer.load_simulation('2D_steady')
air = elmer.load_material('air', sim)
heat_solver = elmer.load_solver('HeatSolver', sim)
output_solver = elmer.load_solver('ResultOutputSolver', sim)
eqn = elmer.Equation(sim, 'main', [heat_solver])
t0 = elmer.InitialCondition(sim, 'T0', {'Temperature': 273.15})

rect = elmer.Body(sim, 'rect', [ph_rect])
rect.material = air
rect.initial_condition = t0
rect.equation = eqn

left = elmer.Boundary(sim, 'left', [bndry_left])
left.fixed_temperature = 293.15
right = elmer.Boundary(sim, 'right', [bndry_right])
right.fixed_temperature = 273.15

sim.write_sif('./simdata')

# run elmer solver
execute.run_elmer_solver('./simdata')
