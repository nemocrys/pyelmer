"""
This example shows the setup of a slightly more complex simulation from
crystal growth.

It investigates the 2d axi-symmetric heat transfer including idealized
radiation and convective cooling modeled with a heat transfer
coefficient.

More information about the case may be found here:
https://www.researchgate.net/publication/340234026_Growth_of_a_tin_crystal
"""

import os
from pyelmer import elmer
from pyelmer import execute
from pyelmer.post import scan_logfile
from objectgmsh import (
    Model,
    Shape,
    MeshControlConstant,
    MeshControlExponential,
    cut,
    factory,
)


###############
# set up working directory
sim_dir = "./simdata"

if not os.path.exists(sim_dir):
    os.mkdir(sim_dir)

###############
# dimensions
cruc_r = 0.06  # crucible radius
cruc_h = 0.03  # crucible height
cruc_hi = 0.015  # crucible height inside
melt_r = 0.025  # melt radius
melt_h = 0.01  # melt height
crys_r = 0.005  # crystal radius
crys_h = 0.1  # crystal height

# geometry modeling using pyelmer.gmsh
model = Model()

# main bodies
crystal = Shape(model, 2, "crystal", [factory.addRectangle(0, 0, 0, crys_r, crys_h)])
melt = Shape(model, 2, "melt", [factory.addRectangle(0, -melt_h, 0, melt_r, melt_h)])

crucible = factory.addRectangle(0, -melt_h - (cruc_h - cruc_hi), 0, cruc_r, cruc_h)
crucible_hole = factory.addRectangle(0, -melt_h, 0, melt_r, cruc_hi)
cut([(2, crucible)], [(2, crucible_hole)])
crucible = Shape(model, 2, "crucible", [crucible])

# create connection between the shapes
crystal.set_interface(melt)
melt.set_interface(crucible)

# detect boundaries
bnd_crystal_out = Shape(
    model, 1, "bnd_crystal_out", [crystal.top_boundary, crystal.right_boundary]
)
bnd_melt = Shape(
    model, 1, "bnd_melt_surf", melt.get_boundaries_in_box([crys_r, melt_r], [0, 0])
)
surfs = [
    crucible.get_boundaries_in_box(
        [melt_r, melt_r], [0, cruc_hi - melt_h], one_only=True
    ),
    crucible.top_boundary,
    crucible.right_boundary,
]
bnd_crucible_outside = Shape(model, 1, "bnd_crucible_outside", surfs)
bnd_crucible_bottom = Shape(model, 1, "bnd_crucible_bottom", [crucible.bottom_boundary])

if_crystal_melt = Shape(model, 1, "if_crystal_melt", crystal.get_interface(melt))

# add physical groups
model.make_physical()

# set mesh constraints
model.deactivate_characteristic_length()
MeshControlConstant(model, 0.005, [crucible, melt])
MeshControlConstant(model, 0.0025, [crystal])
MeshControlExponential(
    model, if_crystal_melt, 0.001, exp=1.7, fact=3, shapes=[crystal, melt, crucible]
)

# create mesh, show, export
model.generate_mesh()
model.show()
model.write_msh(sim_dir + "/case.msh")

###############
# elmer setup

elmer.data_dir = "./data"

sim = elmer.load_simulation("axi-symmetric_steady")

# materials
tin_solid = elmer.load_material("tin_solid", sim)
tin_solid.data.update({"Emissivity": 0.5})
tin_liquid = elmer.load_material("tin_liquid", sim)
tin_liquid.data.update({"Emissivity": 0.5, "Heat Conductivity": 32.0})
graphite = elmer.load_material("graphite_CZ3-R6300", sim)
graphite.data.update({"Emissivity": 0.5, "Heat Conductivity": 236})

# solver, equation
solver_heat = elmer.load_solver("HeatSolver", sim)
solver_output = elmer.load_solver("ResultOutputSolver", sim)
eqn = elmer.Equation(sim, "main", [solver_heat])

# bodies
bdy_crc = elmer.Body(sim, "crucible", [crucible.ph_id])
bdy_crc.material = graphite
bdy_crc.equation = eqn

bdy_melt = elmer.Body(sim, "melt", [melt.ph_id])
bdy_melt.material = tin_liquid
bdy_melt.equation = eqn

bdy_crys = elmer.Body(sim, "crystal", [crystal.ph_id])
bdy_crys.material = tin_solid
bdy_crys.equation = eqn

# boundaries
bndry_bottom = elmer.Boundary(sim, "bottom", [bnd_crucible_bottom.ph_id])
bndry_bottom.data.update({"Heat Flux BC": True, "Heat Flux": 5314})

bndry_crys_melt = elmer.Boundary(sim, "crys_melt", [if_crystal_melt.ph_id])
bndry_crys_melt.data.update({"Heat Flux BC": True, "Heat Flux": 21289})

bndry_surfaces = elmer.Boundary(
    sim, "surfaces", [bnd_crystal_out.ph_id, bnd_melt.ph_id, bnd_crucible_outside.ph_id]
)
bndry_surfaces.data.update({"External Temperature": 300.0})
bndry_surfaces.data.update({"Heat Transfer Coefficient": 3.5})
bndry_surfaces.data.update({"Radiation": "Idealized"})

# alternatively you could load the boundaries from file:

# bndry_bottom = elmer.load_boundary("fixed_heatflux_bottom", sim)
# bndry_bottom.geo_ids = [bnd_crucible_bottom.ph_id]
# bndry_crys_melt = elmer.load_boundary("fixed_heatflux_crys_melt", sim)
# bndry_crys_melt = [if_crystal_melt.ph_id]
# bndry_surfaces = elmer.load_boundary("heat_transfer_radiation_idealized", sim)
# bndry_surfaces.geo_ids = [bnd_crystal_out.ph_id, bnd_melt.ph_id, bnd_crucible_outside.ph_id]

# export
sim.write_startinfo(sim_dir)
sim.write_sif(sim_dir)

##############
# execute ElmerGrid & ElmerSolver
execute.run_elmer_grid(sim_dir, "case.msh")
execute.run_elmer_solver(sim_dir)

###############
# scan log for errors and warnings
err, warn, stats = scan_logfile(sim_dir)
print("Errors:", err)
print("Warnings:", warn)
print("Statistics:", stats)
