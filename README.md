# PyElmer

## Project description

The pyelmer package provides a simple object-oriented way to set up [Elmer FEM](http://www.elmerfem.org/) simulations from python.

Some utility-functions for pre-processing with [pygmsh](https://pypi.org/project/pygmsh/), execution of ElmerGrid and ElmerSolver, and some post-processing routines are provided.

## Installation

The package is not yet listed on [pypi](https://pypi.org/). To install pyelmer you need to clone this repository and run the command
```
pip install -e ./
```
in the project directory.

## Basic principles

The basic working principle of pyelmer is the representation of the entries in the sif-file using dictionaries. Each section of the sif-file is represented by instances of the classes

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

# write sif-file
sim.write_sif('./simulation_directory/')
```

## Examples

### Set up simulation from scratch

To set up a simulation (using an existing geometry) run
```python
import pyelmer.sif as elmer

simulation = elmer.Simulation()
simulation.settings = {
    'Coordinate System': 'Axi Symmetric',
    'Simulation Type': 'Steady state'
}

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

# body with id 1 in mesh
atm = elmer.Body(simulation, 'atmosphere', [1])
# alternatively
atm = elmer.Body(simulation, 'atmosphere')
atm.body_ids = [1]
atm.material = air

# body with id 2, 3 in mesh
wtr = elmer.Body(simulation, 'body2', [2, 3])
wtr.material = water

# To be continued...
```

## Data
TODO