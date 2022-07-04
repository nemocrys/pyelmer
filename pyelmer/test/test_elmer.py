import pytest
import yaml
import gmsh
import os
import re
from objectgmsh import add_physical_group, get_boundaries_in_box
from pyelmer import elmer

elmer.data_dir = "./test_data"

###############
# set up working directory
file_dir = os.path.dirname(__file__)
sim_dir = file_dir + "/test_simdata"

if not os.path.exists(sim_dir):
    os.mkdir(sim_dir)


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
    test_initial_condtion = elmer.InitialCondition(sim, "test_initial_condition")
    test_solver = elmer.Solver(sim, "test_solver")
    test_eqn = elmer.Equation(sim, "main", [test_solver])

    # setup body object
    test_body = elmer.Body(sim, "test_body", [ph_test_body])

    # initialize body data
    test_body.material = elmer.Material(sim, "test_material")
    test_body.initial_condition = test_initial_condtion
    test_body.equation = test_eqn

    assert test_body.get_data() == {
        "Target Bodies(1)": "1",
        "Equation": "0  ! main",
        "Initial Condition": "0  ! test_initial_condition",
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
    test_boundary.data.update({"test_parameter": 1.0})
    assert test_boundary.get_data() == {
        "Target Boundaries(1)": "1",
        "test_parameter": 1.0,
    }

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
    data_test_material = {"test_parameter": 1.0}
    test_material = elmer.Material(sim, "test_material", data=data_test_material)
    assert test_material.get_data() == {"test_parameter": 1.0}
    # initialize empty test material
    empty_test_material = elmer.Material(sim, "empty_test_material")
    assert empty_test_material.get_data() == {}


def test_body_force():
    # setup simulation
    sim = elmer.Simulation()
    # initialize data for test body force
    data_test_body_force = {"test_parameter": 1.0}
    test_material = elmer.BodyForce(sim, "test_body_force", data=data_test_body_force)
    assert test_material.get_data() == {"test_parameter": 1.0}
    # initialize empty test body force
    empty_test_material = elmer.BodyForce(sim, "empty_test_body_force")
    assert empty_test_material.get_data() == {}


def test_initial_condition():
    # setup simulation
    sim = elmer.Simulation()
    # initialize test initial condition
    test_initial_condition = elmer.InitialCondition(
        sim, "test_initial_condition", {"test_parameter": 1.0}
    )
    assert test_initial_condition.get_data() == {"test_parameter": 1.0}
    # initialize empty test initial condition
    empty_test_initial_condition = elmer.InitialCondition(
        sim, "empty_test_initial_condition"
    )
    assert empty_test_initial_condition.get_data() == {}


def test_solver():
    # setup simulation
    sim = elmer.Simulation()

    # initialize test solver
    data_test_solver = {"test_parameter": 1.0}
    test_solver = elmer.Solver(sim, "test_solver", data=data_test_solver)
    assert test_solver.get_data() == {"test_parameter": 1.0}

    # initialize empty test solver
    empty_test_solver = elmer.Solver(sim, "empty_test_solver")
    assert empty_test_solver.get_data() == {}


def test_equation():
    # setup simulation
    sim = elmer.Simulation()

    test_solver = elmer.Solver(sim, "test_solver")

    test_eqn = elmer.Equation(sim, "main", [test_solver])
    assert test_eqn.get_data() == {"Active Solvers(1)": "0   ! test_solver, "}


def test_load_simulation():
    with open(file_dir + "/test_data/simulations.yml") as f:
        settings = yaml.safe_load(f)["test_simulation"]
    assert (
        elmer.load_simulation(
            "test_simulation", file_dir + "/test_data/simulations.yml"
        ).settings
        == settings
    )


def test_load_material():
    with open(file_dir + "/test_data/materials.yml") as f:
        data = yaml.safe_load(f)["test_material"]
    sim = elmer.Simulation()
    assert (
        elmer.load_material(
            "test_material", sim, file_dir + "/test_data/materials.yml"
        ).get_data()
        == data
    )


def test_load_solver():
    with open(file_dir + "/test_data/solvers.yml") as f:
        data = yaml.safe_load(f)["test_solver"]
    sim = elmer.Simulation()
    assert (
        elmer.load_solver(
            "test_solver", sim, file_dir + "/test_data/solvers.yml"
        ).get_data()
        == data
    )


def test_load_body_force():
    with open(file_dir + "/test_data/body_forces.yml") as f:
        data = yaml.safe_load(f)["test_body_force"]
    sim = elmer.Simulation()
    assert (
        elmer.load_body_force(
            "test_body_force", sim, file_dir + "/test_data/body_forces.yml"
        ).get_data()
        == data
    )


def test_load_boundary():
    with open(file_dir + "/test_data/boundaries.yml") as f:
        data = yaml.safe_load(f)["test_boundary"]
    sim = elmer.Simulation()
    assert elmer.load_boundary(
        "test_boundary", sim, file_dir + "/test_data/boundaries.yml"
    ).get_data() == {
        "Target Boundaries(0)": "",
        "test_parameter": data["test_parameter"],
    }


def test_load_initial_condition():
    with open(file_dir + "/test_data/initial_conditions.yml") as f:
        data = yaml.safe_load(f)["test_initial_condition"]
    sim = elmer.Simulation()
    assert (
        elmer.load_initial_condition(
            "test_initial_condition",
            sim,
            file_dir + "/test_data/initial_conditions.yml",
        ).get_data()
        == data
    )


def test_write_sif():
    # setup gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("test-2d")
    factory = gmsh.model.occ

    # create test domain
    test_domain = factory.addRectangle(0, 0, 0, 1, 1)
    factory.synchronize()
    ph_test_body = add_physical_group(2, [test_domain], "test_material")

    # setup simulation
    sim = elmer.load_simulation(
        "test_simulation", file_dir + "/test_data/simulations.yml"
    )

    test_solver = elmer.load_solver(
        "test_solver", sim, file_dir + "/test_data/solvers.yml"
    )

    test_post_solver = elmer.load_solver(
        "test_post_solver", sim, file_dir + "/test_data/solvers.yml"
    )

    test_eqn = elmer.Equation(sim, "test_equation", [test_solver])

    test_initial_condtion = elmer.load_initial_condition(
        "test_initial_condition", sim, file_dir + "/test_data/initial_conditions.yml"
    )

    test_body_force = elmer.load_body_force(
        "test_body_force", sim, file_dir + "/test_data/body_forces.yml"
    )
    # setup body object
    test_body = elmer.Body(sim, "test_body", [ph_test_body])

    # initialize body data
    test_material = elmer.load_material(
        "test_material", sim, file_dir + "/test_data/materials.yml"
    )
    test_body.material = test_material
    test_body.initial_condition = test_initial_condtion
    test_body.equation = test_eqn

    # detect boundary
    line = get_boundaries_in_box(0, 0, 0, 1, 0, 0, 2, test_domain)
    ph_test_boundary = add_physical_group(1, [line], "test_boundary")

    # initialize test boundary
    test_boundary = elmer.Boundary(sim, "test_boundary", [ph_test_boundary])
    test_boundary.data.update({"test_parameter": "1.0"})

    sim.write_startinfo(sim_dir)
    sim.write_sif(sim_dir)

    # check format of sif file
    simulation = {}
    constants = {}
    equations = {}
    solvers = {}
    materials = {}
    bodies = {}
    boundaries = {}
    body_forces = {}
    initial_conditions = {}
    names = [
        "Simulation",
        "Constants",
        "Equation",
        "Solver",
        "Material",
        "Body",
        "Boundary Condition",
        "Body Force",
        "Initial Condition",
    ]
    objects = [
        simulation,
        constants,
        equations,
        solvers,
        materials,
        bodies,
        boundaries,
        body_forces,
        initial_conditions,
    ]

    with open(sim_dir + "/case.sif", "r") as f:
        read_object = False
        write_index = None
        write_name = None
        checked_name = None
        object_name = None
        while True:
            line = f.readline()

            if not line:
                break

            line = line.strip().strip("\n")
            if line in names:
                read_object = True
                for i, name in enumerate(names):
                    if name == line.strip("\n"):
                        write_index = i
                continue

            if line == "End":
                read_object = False
                checked_name = None
                continue

            if checked_name is not None and not read_object:
                read_object = True
                if any(map(str.isdigit, line)):
                    object_name, object_number, _ = map(
                        str.strip, re.split("(\d+)", line)
                    )
                    for i, name in enumerate(names):
                        if name == object_name:
                            write_index = i
                            write_name = name
                    objects[write_index].update({checked_name: {}})
                    continue

            if line.startswith("!"):
                checked_name = line.strip("!").strip(" ")
            if read_object:
                key, value = line.strip(" ").split(" = ")
                if write_name == "Equation":
                    value = value + " "
                if checked_name is not None:
                    objects[write_index][checked_name].update({key: value})
                else:
                    objects[write_index].update({key: value})

    assert equations["test_equation"] == test_eqn.get_data()
    assert solvers["test_solver"] == test_solver.get_data()
    assert solvers["test_post_solver"] == test_post_solver.get_data()
    assert materials["test_material"] == test_material.get_data()
    assert bodies["test_body"] == test_body.get_data()
    assert boundaries["test_boundary"] == test_boundary.get_data()
    assert body_forces["test_body_force"] == test_body_force.get_data()
    assert (
        initial_conditions["test_initial_condition"] == test_initial_condtion.get_data()
    )
