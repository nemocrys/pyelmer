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
