import pytest
import yaml
import gmsh
import os
from pyelmer.gmsh import add_physical_group, get_boundaries_in_box
from pyelmer import elmer

elmer.data_dir = "./test_data"

###############
# set up working directory
sim_dir = "./test_simdata"
file_dir = os.path.dirname(__file__)

if not os.path.exists(sim_dir):
    os.mkdir(sim_dir)


def test_simulation():
    # TO DO
    pass


def test_body():
    # setup gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("test-2d")
    factory = gmsh.model.occ

    # create test domain
    test_domain = factory.addRectangle(0, 0, 0, 1, 1)
    ph_test_body = add_physical_group(2, [test_domain], "test_body")

    # setup simulation
    sim = elmer.Simulation()
    test_initial_condtion = elmer.InitialCondition(sim, "T0", {"Temperature": 1.0})
    test_solver = elmer.Solver(
        sim,
        "test_solver",
        data={
            "Equation": "HeatSolver",
            "Procedure": '"HeatSolve" "HeatSolver"',
            "Variable": '"Temperature"',
            "Variable Dofs": 1,
            "Calculate Loads": True,
            "Exec Solver": "Always",
            "Nonlinear System Convergence Tolerance": 1e-06,
            "Nonlinear System Max Iterations": 1000,
            "Nonlinear System Relaxation Factor": 0.7,
            "Steady State Convergence Tolerance": 1e-06,
            "Stabilize": True,
            "Optimize Bandwidth": True,
            "Linear System Solver": "Iterative",
            "Linear System Iterative Method": "BiCGStab",
            "Linear System Max Iterations": 1000,
            "Linear System Preconditioning": "ILU",
            "Linear System Precondition Recompute": 1,
            "Linear System Convergence Tolerance": 1e-08,
            "Linear System Abort Not Converged": True,
            "Linear System Residual Output": 1,
            "Smart Heater Control After Tolerance": 0.0001,
        },
    )
    test_eqn = elmer.Equation(sim, "main", [test_solver])

    # setup body object
    test_body = elmer.Body(sim, "test_body", [ph_test_body])

    # initialize body data
    test_material = {"Density": 1.0, "Heat Capacity": 1.0, "Heat Conductivity": 1.0}
    test_body.material = elmer.Material(sim, "test_material", data=test_material)
    test_body.initial_condition = test_initial_condtion
    test_body.equation = test_eqn

    assert test_body.get_data() == {
        "Target Bodies(1)": "1",
        "Equation": "0  ! main",
        "Initial Condition": "0  ! T0",
        "Material": "0  ! test_material",
    }

    # initialize empty body
    empty_test_body = elmer.Body(sim, "empty_test_body")
    assert empty_test_body.get_data() == {"Target Bodies(0)": ""}


def test_boundary():
    # setup gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("test-2d")
    factory = gmsh.model.occ

    # create test domain
    test_domain = factory.addRectangle(0, 0, 0, 1, 1)
    factory.synchronize()

    # detect boundary
    line = get_boundaries_in_box(0, 0, 0, 1, 0, 0, 2, test_domain)
    ph_test_boundary = add_physical_group(1, [line], "test_boundary")

    # setup simulation
    sim = elmer.Simulation()
    # initialize test boundary
    test_boundary = elmer.Boundary(sim, "test_boundary", [ph_test_boundary])
    test_boundary.data.update({"Temperature": 1.0})
    assert test_boundary.get_data() == {"Target Boundaries(1)": "1", "Temperature": 1.0}

    # initialize flux test boundary
    flux_test_boundary = elmer.Boundary(sim, "flux_test_boundary", [ph_test_boundary])
    flux_test_boundary.data.update({"Heat Flux BC": True, "Heat Flux": 1})
    assert flux_test_boundary.get_data() == {
        "Target Boundaries(1)": "1",
        "Heat Flux BC": True,
        "Heat Flux": 1.0,
    }

    # initialize empty test boundary
    empty_test_boundary = elmer.Boundary(sim, "empty_test_boundary")
    assert empty_test_boundary.get_data() == {"Target Boundaries(0)": ""}


def test_material():
    # setup simulation
    sim = elmer.Simulation()
    # initialize data for test material
    data_test_material = {
        "Density": 1.0,
        "Electric Conductivity": 1.0,
        "Emissivity": 1.0,
        "Heat Capacity": 1.0,
        "Heat Conductivity": 1.0,
        "Relative Permeability": 1,
        "Relative Permittivity": 1,
        "Solid": "Logical True",
        "Melting Point": 1,
        "Latent Heat": 1,
    }
    test_material = elmer.Material(sim, "test_material", data=data_test_material)
    assert test_material.get_data() == {
        "Density": 1.0,
        "Electric Conductivity": 1.0,
        "Emissivity": 1.0,
        "Heat Capacity": 1.0,
        "Heat Conductivity": 1.0,
        "Relative Permeability": 1,
        "Relative Permittivity": 1,
        "Solid": "Logical True",
        "Melting Point": 1,
        "Latent Heat": 1,
    }
    # initialize empty test material
    empty_test_material = elmer.Material(sim, "empty_test_material")
    assert empty_test_material.get_data() == {}


def test_body_force():
    # TO DO
    pass


def test_initial_condition():
    # setup simulation
    sim = elmer.Simulation()
    # initialize test initial condition
    test_initial_condition = elmer.InitialCondition(
        sim, "test_initial_condition", {"Temperature": 1.0}
    )
    assert test_initial_condition.get_data() == {"Temperature": 1.0}
    # initialize empty test initial condition
    empty_test_initial_condition = elmer.InitialCondition(
        sim, "empty_test_initial_condition"
    )
    assert empty_test_initial_condition.get_data() == {}


def test_solver():
    # setup simulation
    sim = elmer.Simulation()

    # initialize test solver
    data_test_solver = {
        "Equation": "HeatSolver",
        "Procedure": '"HeatSolve" "HeatSolver"',
        "Variable": '"Temperature"',
        "Variable Dofs": 1,
        "Calculate Loads": True,
        "Exec Solver": "Always",
        "Nonlinear System Convergence Tolerance": 1e-06,
        "Nonlinear System Max Iterations": 1000,
        "Nonlinear System Relaxation Factor": 0.7,
        "Steady State Convergence Tolerance": 1e-06,
        "Stabilize": True,
        "Optimize Bandwidth": True,
        "Linear System Solver": "Iterative",
        "Linear System Iterative Method": "BiCGStab",
        "Linear System Max Iterations": 1000,
        "Linear System Preconditioning": "ILU",
        "Linear System Precondition Recompute": 1,
        "Linear System Convergence Tolerance": 1e-08,
        "Linear System Abort Not Converged": True,
        "Linear System Residual Output": 1,
        "Smart Heater Control After Tolerance": 0.0001,
    }
    test_solver = elmer.Solver(sim, "test_solver", data=data_test_solver)
    assert test_solver.get_data() == {
        "Equation": "HeatSolver",
        "Procedure": '"HeatSolve" "HeatSolver"',
        "Variable": '"Temperature"',
        "Variable Dofs": 1,
        "Calculate Loads": True,
        "Exec Solver": "Always",
        "Nonlinear System Convergence Tolerance": 1e-06,
        "Nonlinear System Max Iterations": 1000,
        "Nonlinear System Relaxation Factor": 0.7,
        "Steady State Convergence Tolerance": 1e-06,
        "Stabilize": True,
        "Optimize Bandwidth": True,
        "Linear System Solver": "Iterative",
        "Linear System Iterative Method": "BiCGStab",
        "Linear System Max Iterations": 1000,
        "Linear System Preconditioning": "ILU",
        "Linear System Precondition Recompute": 1,
        "Linear System Convergence Tolerance": 1e-08,
        "Linear System Abort Not Converged": True,
        "Linear System Residual Output": 1,
        "Smart Heater Control After Tolerance": 0.0001,
    }

    # initialize empty test solver
    empty_test_solver = elmer.Solver(sim, "empty_test_solver")
    assert empty_test_solver.get_data() == {}


def test_equation():
    # setup simulation
    sim = elmer.Simulation()

    data_test_solver = {
        "Equation": "HeatSolver",
        "Procedure": '"HeatSolve" "HeatSolver"',
        "Variable": '"Temperature"',
        "Variable Dofs": 1,
        "Calculate Loads": True,
        "Exec Solver": "Always",
        "Nonlinear System Convergence Tolerance": 1e-06,
        "Nonlinear System Max Iterations": 1000,
        "Nonlinear System Relaxation Factor": 0.7,
        "Steady State Convergence Tolerance": 1e-06,
        "Stabilize": True,
        "Optimize Bandwidth": True,
        "Linear System Solver": "Iterative",
        "Linear System Iterative Method": "BiCGStab",
        "Linear System Max Iterations": 1000,
        "Linear System Preconditioning": "ILU",
        "Linear System Precondition Recompute": 1,
        "Linear System Convergence Tolerance": 1e-08,
        "Linear System Abort Not Converged": True,
        "Linear System Residual Output": 1,
        "Smart Heater Control After Tolerance": 0.0001,
    }
    test_solver = elmer.Solver(sim, "test_solver", data=data_test_solver)

    test_eqn = elmer.Equation(sim, "main", [test_solver])
    assert test_eqn.get_data() == {"Active Solvers(1)": "0   ! test_solver, "}


@pytest.mark.parametrize(
    "name",
    ["2D_steady", "2D_transient", "axi-symmetric_steady", "axi-symmetric_transient",],
)
def test_load_simulation(name):
    with open(file_dir + "/test_data/simulations.yml") as f:
        settings = yaml.safe_load(f)[name]
    assert (
        elmer.load_simulation(name, file_dir + "/test_data/simulations.yml").settings
        == settings
    )


@pytest.mark.parametrize("material", ["air", "water", "tin_liquid", "tin_solid",])
def test_load_material(material):
    with open(file_dir + "/test_data/materials.yml") as f:
        data = yaml.safe_load(f)[material]
    sim = elmer.Simulation()
    assert (
        elmer.load_material(
            material, sim, file_dir + "/test_data/materials.yml"
        ).get_data()
        == data
    )


@pytest.mark.parametrize(
    "solver",
    [
        "ThermoElectricSolver",
        "HeatSolver",
        "MagnetoDynamics2DHarmonic",
        "MagnetoDynamicsCalcFields",
        "StatMagSolver",
        "SaveMaterials",
        "ResultOutputSolver",
        "FluxSolver",
        "SaveScalars",
        "SaveLine",
        "SteadyPhaseChange",
    ],
)
def test_load_solver(solver):
    with open(file_dir + "/test_data/solvers.yml") as f:
        data = yaml.safe_load(f)[solver]
    sim = elmer.Simulation()
    assert (
        elmer.load_solver(solver, sim, file_dir + "/test_data/solvers.yml").get_data()
        == data
    )
