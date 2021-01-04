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
import gmsh
from pyelmer import elmer
from pyelmer import execute
from pyelmer.post import scan_logfile
from pyelmer.gmsh_utils import add_physical_group, get_boundaries_in_box


###############
# set up working directory
sim_dir = './simdata'

if not os.path.exists(sim_dir):
    os.mkdir(sim_dir)

###############
# geometry modeling using gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add('heat-transfer-2d')
factory = gmsh.model.occ

# dimensions
cruc_r = 0.06
cruc_h = 0.03
cruc_hd = 0.015
melt_r = 0.025
melt_h = 0.01
crys_r = 0.005
crys_h = 0.1

# main bodies
crys = factory.addRectangle(0, 0, 0, crys_r, crys_h)
melt = factory.addRectangle(0, -melt_h, 0, melt_r, melt_h)
crucin = factory.addRectangle(0, -melt_h-(cruc_h-cruc_hd), 0, melt_r, cruc_h-cruc_hd)
crucout = factory.addRectangle(melt_r, -melt_h-(cruc_h-cruc_hd), 0, cruc_r-melt_r, cruc_h)

# create connection between the two bodies
factory.synchronize()
factory.fragment([(2, melt)], [(2, crys), (2, crucin), (2, crucout)])
factory.fragment([(2, crucin)], [(2, crucout)])

# add physical groups
factory.synchronize()
ph_crys = add_physical_group(2, [crys], 'crys')
ph_melt = add_physical_group(2, [melt], 'melt')
ph_cruc = add_physical_group(2, [crucin, crucout], 'cruc')

# detect boundaries 
crucbotin = get_boundaries_in_box(0, -(cruc_hd+melt_h), 0, cruc_r, -(cruc_hd+melt_h), 0, 2, crucin)
crucbotout = get_boundaries_in_box(0, -(cruc_hd+melt_h), 0, cruc_r, -(cruc_hd+melt_h), 0, 2, crucout)
crucside = get_boundaries_in_box(cruc_r, -(cruc_hd+melt_h), 0, cruc_r, -(cruc_hd+melt_h)+cruc_h, 0, 2, crucout)
cructop = get_boundaries_in_box(melt_r, -(cruc_hd+melt_h)+cruc_h, 0, cruc_r, -(cruc_hd+melt_h)+cruc_h, 0, 2, crucout)
crucin = get_boundaries_in_box(melt_r, 0, 0, melt_r, -(cruc_hd+melt_h)+cruc_h, 0, 2, crucout)
melttop = get_boundaries_in_box(crys_r, 0, 0, melt_r, 0, 0, 2, melt)
crysint = get_boundaries_in_box(0, 0, 0, crys_r, 0, 0, 2, crys)
cryssurf = get_boundaries_in_box(crys_r, 0, 0, crys_r, crys_h, 0, 2, crys)
crystop = get_boundaries_in_box(0, crys_h, 0, crys_r, crys_h, 0, 2, crys)

ph_bottom = add_physical_group(1, [crucbotin, crucbotout], 'bottom')
ph_rad = add_physical_group(1, [crucside, cructop, crucin, melttop, cryssurf, crystop], 'rad')
ph_crysint = add_physical_group(1, [crysint], 'crysint')

# create mesh
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.002)
gmsh.model.mesh.generate(2)

# show mesh & export
gmsh.fltk.run()
gmsh.write(sim_dir + '/case.msh2')

###############
# elmer setup
sim = elmer.load_simulation('axi-symmetric_steady')

# materials
tin_solid = elmer.load_material('tin_solid', sim)
tin_solid.data.update({'Emissivity': 0.5})
tin_liquid = elmer.load_material('tin_liquid', sim)
tin_liquid.data.update({'Emissivity': 0.5,
                        'Heat Conductivity': 32.})
graphite = elmer.load_material('graphite_CZ3-R6300', sim)
graphite.data.update({'Emissivity': 0.5,
                      'Heat Conductivity': 236})

# solver, equation
solver_heat = elmer.load_solver('HeatSolver', sim)
solver_output = elmer.load_solver('ResultOutputSolver', sim)
eqn = elmer.Equation(sim, 'main', [solver_heat])

# initial condition
T0 = elmer.InitialCondition(sim, 'T0', {'Temperature': 273.15})

# bodies
bdy_crc = elmer.Body(sim, 'crucible', [ph_cruc])
bdy_crc.material = graphite
bdy_crc.initial_condition = T0
bdy_crc.equation = eqn

bdy_melt = elmer.Body(sim, 'melt', [ph_melt])
bdy_melt.material = tin_liquid
bdy_melt.initial_condition = T0
bdy_melt.equation = eqn

bdy_crys = elmer.Body(sim, 'crystal', [ph_crys])
bdy_crys.material = tin_solid
bdy_crys.initial_condition = T0
bdy_crys.equation = eqn

# boundaries
bndry_bottom = elmer.Boundary(sim, 'bottom', [ph_bottom])
bndry_bottom.data.update({'Heat Flux BC': True,
                          'Heat Flux': 'real 5314'})

bndry_crysint = elmer.Boundary(sim, 'crysint', [ph_crysint])
bndry_crysint.data.update({'Heat Flux BC': True,
                           'Heat Flux': 'real 21289'})

bndry_rad = elmer.Boundary(sim, 'rad', [ph_rad])
bndry_rad.heat_transfer_coefficient = 3.5
bndry_rad.T_ext = 300.0
bndry_rad.data.update({'Radiation': 'Idealized'})

# export
sim.write_startinfo(sim_dir)
sim.write_sif(sim_dir)

##############
# execute ElmerGrid & ElmerSolver
execute.run_elmer_grid(sim_dir, 'case.msh2')
execute.run_elmer_solver(sim_dir)

###############
# scan log for errors and warnings
err, warn, stats = scan_logfile(sim_dir)
print('Errors:', err)
print('Warnings:', warn)
print('Statistics:', stats)
