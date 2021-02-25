"""
This is an alternative version of the example heat_transfer_2d.py

Instead of loading pre defined solvers and material data, the whole case
is set up manually.
"""
import os
import gmsh
from pyelmer import elmer
from pyelmer import execute
from pyelmer.post import scan_logfile
from pyelmer.gmsh_utils import add_physical_group, get_boundaries_in_box


###############
# set up working directory
sim_dir = "./simdata"

if not os.path.exists(sim_dir):
    os.mkdir(sim_dir)

###############
# geometry modeling using gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("heat-transfer-2d")
factory = gmsh.model.occ

# main bodies
water = factory.addRectangle(0, 0, 0, 1, 1)
air = factory.addRectangle(0, 1, 0, 1, 1)

# create connection between the two bodies
factory.synchronize()
factory.fragment([(2, water)], [(2, air)])

# add physical groups
factory.synchronize()
ph_water = add_physical_group(2, [water], "water")
ph_air = add_physical_group(2, [air], "air")

# detect boundaries
line = get_boundaries_in_box(0, 0, 0, 1, 0, 0, 2, water)
ph_bottom = add_physical_group(1, [line], "bottom")
line = get_boundaries_in_box(0, 2, 0, 1, 2, 0, 2, air)
ph_top = add_physical_group(1, [line], "top")

# create mesh
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.1)
gmsh.model.mesh.generate(2)

# show mesh & export
gmsh.fltk.run()  # comment this line out if your system doesn't support the gmsh GUI
gmsh.write(sim_dir + "/case.msh2")  # use legacy file format msh2 for elmer grid

###############
# elmer setup
sim = elmer.Simulation()
sim.settings = {"Coordinate System": "Cartesian 2D", "Simulation Type": "Steady state"}

air = elmer.Material(sim, "air")
air.data = {"Density": 1.1885, "Heat Capacity": 1006.4, "Heat Conductivity": 0.025873}
water = elmer.Material(sim, "water")
water.data = {"Density": 1000.0, "Heat Capacity": 4182.0, "Heat Conductivity": 0.6}

solver_heat = elmer.Solver(sim, "heat_solver")
solver_heat.data = {
    "Equation": "HeatSolver",
    "Procedure": '"HeatSolve" "HeatSolver"',
    "Variable": '"Temperature"',
    "Variable Dofs": 1,
}
solver_output = elmer.Solver(sim, "output_solver")
solver_output.data = {
    "Exec Solver": "After timestep",
    "Equation": "ResultOutputSolver",
    "Procedure": '"ResultOutputSolve" "ResultOutputSolver"',
}
eqn = elmer.Equation(sim, "main", [solver_heat])

T0 = elmer.InitialCondition(sim, "T0", {"Temperature": 273.15})

bdy_water = elmer.Body(sim, "water", [ph_water])
bdy_water.material = water
bdy_water.initial_condition = T0
bdy_water.equation = eqn

bdy_air = elmer.Body(sim, "air", [ph_air])
bdy_air.material = air
bdy_air.initial_condition = T0
bdy_air.equation = eqn

bndry_bottom = elmer.Boundary(sim, "bottom", [ph_bottom])
bndry_bottom.fixed_temperature = 353.15  # 80 °C
bndry_top = elmer.Boundary(sim, "top", [ph_top])
bndry_top.fixed_temperature = 293.15  # 20 °C

sim.write_startinfo(sim_dir)
sim.write_sif(sim_dir)

##############
# execute ElmerGrid & ElmerSolver
execute.run_elmer_grid(sim_dir, "case.msh2")
execute.run_elmer_solver(sim_dir)

###############
# scan log for errors and warnings
err, warn, stats = scan_logfile(sim_dir)
print("Errors:", err)
print("Warnings:", warn)
print("Statistics:", stats)
