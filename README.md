# PyElmer

## Project description

The pyelmer package provides a simple object-oriented way to set up [Elmer FEM](http://www.elmerfem.org/) simulations from python.

Some utility-functions for pre-processing using the [gmsh python API](https://pypi.org/project/gmsh/), execution of ElmerGrid and ElmerSolver, and some post-processing routines are provided. Some default simulation settings, solvers, and materials are available.

## Prerequisites

Pyelmer is written in Python 3.7. To run simulations, an Elmer executable is required. As pyelmer was developed to be used with gmsh, an installation of this package is required. Simulation settings, solver, and materials are stored in yaml-files. Therefore pyelmer depends on pyyaml.

```
pip install --upgrade gmsh
pip install --upgrade pyyaml
```

## Installation

The package is not yet listed on [pypi](https://pypi.org/), so pip install pyelmer does not work. To install pyelmer you need to clone this repository and install it in development mode. This is done with the following commands:

```
git clone https://github.com/nemocrys/pyelmer
cd pyelmer
pip install -e ./
```

Alternatively, you can build your own wheel and install it with the tool of your choice (more information about this can be found [here](https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html)).

## Basic principles

The basic working principle of pyelmer is the representation of sif-file entries in dictionaries. Each section of the sif-file is represented by instances of the classes

- *Solver*
- *Equation*
- *Material*
- *Body*
- *Boundary*
- *BodyForce*
- *InitialCondition*

The parameters of e.g. a material are stored in

```python
material.data = {
    'Density': 1.1885,
    'Heat Capacity': 1006.4,
    'Heat Conductivity': 0.025873
}
```

An object of the class *Simulation* is used to manage all members of the sif-file:

```python
import pyelmer.sif as elmer

# simulation object
sim = elmer.Simulation()

# create material and add it to sim
air = elmer.Material(sim, 'air')
air.data = {
    # material data here
}

# add solvers, equations, bodies, ...
heat_solver = elmer.Solver(sim, 'heat')
heat_solver.data = {
    # solver data here
}
# ...

# write sif-file
sim.write_sif('./simulation_directory/')
```

## Examples

### Set up simulation from scratch

```python
from pyelmer import elmer

# set up simulation
sim = elmer.Simulation()
sim.settings = {
    'Coordinate System': 'Axi Symmetric',
    'Simulation Type': 'Steady state'
}
# materials
air = elmer.Material(sim, 'air')
air.data = {
    'Density': 1.1885,
    'Heat Capacity': 1006.4,
    'Heat Conductivity': 0.025873
}
water = elmer.Material(sim, 'water')
water.data = {
    'Density': 1000,
    'Heat Capacity': 4184,
    'Heat Conductivity': 0.6
}
# solvers
heat_solver = elmer.Solver(sim, 'heat')
heat_solver.data = {
    'Equation': 'Heat Equation',
    'Procedure': '"HeatSolve" "HeatSolver"',
    'Variable': '"Temperature"',
    'Variable Dofs': 1
}
output_solver = elmer.Solver(sim, 'output')
output_solver.data = {
    'Exec Solver': 'After timestep',
    'Equation': '"ResultOutput"',
    'Procedure': '"ResultOutputSolve" "ResultOutputSolver"',
    'VTU Format': True
}
# equation
main_eqn = elmer.Equation(sim, 'main', [heat_solver, output_solver])

# initial condition
t0 = elmer.InitialCondition(sim, 'T0')
t0.data = {
    'Temperature': 293.15
}

# body force
heating = elmer.BodyForce(sim, 'heating')
heating.heat_source = 100  # W / kg
# or alternatively:
heating.data = {
    'Heat Source': 100
}

# body with id 1 in mesh
atm = elmer.Body(sim, 'atmosphere')
atm.body_ids = [1]
# or alternatively:
atm = elmer.Body(sim, 'atmosphere', [1])
# set dependencies
atm.material = air
atm.equation = main_eqn
atm.initial_condition = t0

# body with id 2, 3 in mesh
wtr = elmer.Body(sim, 'body2', [2, 3])
wtr.material = water
wtr.equation = main_eqn
wtr.initial_condition = t0
wtr.body_force = heating

# boundary with id 1, 2 in mesh
t1 = elmer.Boundary(sim, 'T1', [1, 2])
t1.fixed_temperature = 293.15
# or:
t1.data ={
    'Temperature': 293.15
}

# write sif file
sim.write_sif('./simulation_directory/')
```

### Execute simulation

```python
from pyelmer import execute

execute.run_elmer_solver('./simulation_directory/', 'ElmerSolver.exe')
```

### Load settings

It is possible to load some default setups for simulations, solvers and materials:

```python
import pyelmer.sif as elmer

# load settings for steady state axi symmetric simulation
sim = elmer.load_simulation('axi-symmetric_steady')
print(sim.settings)

# load heat solver and add it to sim
heat_solver = elmer.load_solver('HeatSolver', sim)
print(heat_solver.data)

# load material and add it to sim
air = elmer.load_material('air', sim)
print(air.data)
```

Unfortunately, there are only very few setups available. You can find them in pyelmer/data.

### Working with gmsh

Pyelmer was developed to facilitate the workflow for geometry generation with gmsh and simulation with elmer. Some utility functions for gmsh (mostly for axi-symmetric cases) are provided.

```python
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

```

## License

Pyelmer is published under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.html).

## Acknowledgements

This package was developed in the [Nemocrys project](https://www.researchgate.net/project/NEMOCRYS-Next-Generation-Multiphysical-Models-for-Crystal-Growth-Processes) which is founded by the [ERC](https://erc.europa.eu/).

## Contribution

Any help to improve this package is very welcome!
