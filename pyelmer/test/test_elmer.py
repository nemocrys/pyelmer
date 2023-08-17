import yaml
import filecmp
import gmsh
import pytest
import os
import re
import tempfile
from objectgmsh import add_physical_group, get_boundaries_in_box
from pyelmer import elmer

elmer.data_dir = os.path.join(os.getcwd(), "test_data")

###############
# set up working directory
file_dir = os.path.dirname(__file__)
sim_dir = os.path.join(file_dir, "test_simdata")

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


def test_component():
    # setup simulation
    sim = elmer.Simulation()
    body1 = elmer.Body(sim, "body1")
    body2 = elmer.Body(sim, "body2")
    boundary1 = elmer.Boundary(sim, "boundary1")
    boundary2 = elmer.Boundary(sim, "boundary2")
    # set id manually for testing purposes, don't do that in a productive environment
    body1.id = 1
    body2.id = 2
    boundary1.id = 3
    boundary2.id = 4

    # initialize data for test component
    test_component = elmer.Component(
        sim,
        "test_component",
        master_bodies=[body1, body2],
        master_boundaries=[boundary1, boundary2],
        data={"test_parameter": 1.0},
    )
    assert test_component.get_data() == {
        "Name": "test_component",
        "Master Bodies(2)": "1 2  ! body1, body2,",
        "Master Boundaries(2)": "3 4  ! boundary1, boundary2,",
        "test_parameter": 1.0,
    }
    empty_test_component = elmer.Component(sim, "empty_test_component")
    assert empty_test_component.get_data() == {
        "Name": "empty_test_component",
    }


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
    with open(os.path.join(file_dir, "test_data", "simulations.yml")) as f:
        settings = yaml.safe_load(f)["test_simulation"]
    assert (
        elmer.load_simulation(
            "test_simulation", os.path.join(file_dir, "test_data", "simulations.yml")
        ).settings
        == settings
    )


def test_load_material():
    with open(os.path.join(file_dir, "test_data", "materials.yml")) as f:
        data = yaml.safe_load(f)["test_material"]
    sim = elmer.Simulation()
    assert (
        elmer.load_material(
            "test_material", sim, os.path.join(file_dir, "test_data", "materials.yml")
        ).get_data()
        == data
    )


def test_load_solver():
    with open(os.path.join(file_dir, "test_data", "solvers.yml")) as f:
        data = yaml.safe_load(f)["test_solver"]
    sim = elmer.Simulation()
    assert (
        elmer.load_solver(
            "test_solver", sim, os.path.join(file_dir, "test_data", "solvers.yml")
        ).get_data()
        == data
    )


def test_load_body_force():
    with open(os.path.join(file_dir, "test_data", "body_forces.yml")) as f:
        data = yaml.safe_load(f)["test_body_force"]
    sim = elmer.Simulation()
    assert (
        elmer.load_body_force(
            "test_body_force",
            sim,
            os.path.join(file_dir, "test_data", "body_forces.yml"),
        ).get_data()
        == data
    )


def test_load_boundary():
    with open(os.path.join(file_dir, "test_data", "boundaries.yml")) as f:
        data = yaml.safe_load(f)["test_boundary"]
    sim = elmer.Simulation()
    assert elmer.load_boundary(
        "test_boundary", sim, os.path.join(file_dir, "test_data", "boundaries.yml")
    ).get_data() == {
        "Target Boundaries(0)": "",
        "test_parameter": data["test_parameter"],
    }


def test_load_initial_condition():
    with open(os.path.join(file_dir, "test_data", "initial_conditions.yml")) as f:
        data = yaml.safe_load(f)["test_initial_condition"]
    sim = elmer.Simulation()
    assert (
        elmer.load_initial_condition(
            "test_initial_condition",
            sim,
            os.path.join(file_dir, "test_data", "initial_conditions.yml"),
        ).get_data()
        == data
    )


def test_duplicate_body():
    sim = elmer.Simulation()

    obj1 = elmer.Body(sim, "duplicate", body_ids=[1], data={"test": 7})
    obj2 = elmer.Body(sim, "duplicate", body_ids=[1], data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.Body(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Body(sim, "duplicate", body_ids=[1])

    with pytest.raises(ValueError):
        elmer.Body(sim, "duplicate", data={"test": 7})


def test_duplicate_boundary():
    sim = elmer.Simulation()

    obj1 = elmer.Boundary(sim, "duplicate", geo_ids=[1], data={"test": 7})
    obj2 = elmer.Boundary(sim, "duplicate", geo_ids=[1], data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.Boundary(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Boundary(sim, "duplicate", geo_ids=[1])

    with pytest.raises(ValueError):
        elmer.Boundary(sim, "duplicate", data={"test": 7})


def test_duplicate_material():
    sim = elmer.Simulation()

    obj1 = elmer.Material(sim, "duplicate", data={"test": 7})
    obj2 = elmer.Material(sim, "duplicate", data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.Material(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Material(sim, "duplicate")

    with pytest.raises(ValueError):
        elmer.Material(sim, "duplicate", data={"test": 8})


def test_duplicate_body_force():
    sim = elmer.Simulation()

    obj1 = elmer.BodyForce(sim, "duplicate", data={"test": 7})
    obj2 = elmer.BodyForce(sim, "duplicate", data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.BodyForce(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.BodyForce(sim, "duplicate")

    with pytest.raises(ValueError):
        elmer.BodyForce(sim, "duplicate", data={"test": 8})


def test_duplicate_initial_condition():
    sim = elmer.Simulation()

    obj1 = elmer.InitialCondition(sim, "duplicate", data={"test": 7})
    obj2 = elmer.InitialCondition(sim, "duplicate", data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.InitialCondition(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.InitialCondition(sim, "duplicate")

    with pytest.raises(ValueError):
        elmer.InitialCondition(sim, "duplicate", data={"test": 8})


def test_duplicate_solver():
    sim = elmer.Simulation()

    obj1 = elmer.Solver(sim, "duplicate", data={"test": 7})
    obj2 = elmer.Solver(sim, "duplicate", data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.Solver(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Solver(sim, "duplicate")

    with pytest.raises(ValueError):
        elmer.Solver(sim, "duplicate", data={"test": 8})


def test_duplicate_equation():
    sim = elmer.Simulation()
    solver1 = elmer.Solver(sim, "solver1")
    solver2 = elmer.Solver(sim, "solver2")

    obj1 = elmer.Equation(sim, "duplicate", [solver1], data={"test": 7})
    obj2 = elmer.Equation(sim, "duplicate", [solver1], data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.Equation(sim, "different", [solver1])
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Equation(sim, "duplicate", [solver1])

    with pytest.raises(ValueError):
        elmer.Equation(sim, "duplicate", [solver2], data={"test": 7})


def test_duplicate_component():
    sim = elmer.Simulation()

    # test with master body
    body = elmer.Body(sim, "test_body")

    obj1 = elmer.Component(sim, "duplicate", master_bodies=[body], data={"test": 7})
    obj2 = elmer.Component(sim, "duplicate", master_bodies=[body], data={"test": 7})
    assert obj1 is obj2

    obj3 = elmer.Component(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Component(sim, "duplicate", master_bodies=[body])

    with pytest.raises(ValueError):
        elmer.Component(sim, "duplicate", data={"test": 7})

    # test with master boundary
    boundary = elmer.Body(sim, "test_boundary")

    obj1 = elmer.Component(
        sim, "duplicate2", master_boundaries=[boundary], data={"test": 7}
    )
    obj2 = elmer.Component(
        sim, "duplicate2", master_boundaries=[boundary], data={"test": 7}
    )
    assert obj1 is obj2

    obj3 = elmer.Component(sim, "different")
    assert obj1 is not obj3

    with pytest.raises(ValueError):
        elmer.Component(sim, "duplicate2", master_boundaries=[boundary])

    with pytest.raises(ValueError):
        elmer.Component(sim, "duplicate2", data={"test": 7})


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
        "test_simulation", os.path.join(file_dir, "test_data", "simulations.yml")
    )

    test_solver = elmer.load_solver(
        "test_solver", sim, os.path.join(file_dir, "test_data", "solvers.yml")
    )

    test_post_solver = elmer.load_solver(
        "test_post_solver", sim, os.path.join(file_dir, "test_data", "solvers.yml")
    )

    test_eqn = elmer.Equation(sim, "test_equation", [test_solver])

    test_initial_condtion = elmer.load_initial_condition(
        "test_initial_condition",
        sim,
        os.path.join(file_dir, "test_data", "initial_conditions.yml"),
    )

    test_body_force = elmer.load_body_force(
        "test_body_force", sim, os.path.join(file_dir, "test_data", "body_forces.yml")
    )
    # setup body object
    test_body = elmer.Body(sim, "test_body", [ph_test_body])

    # initialize body data
    test_material = elmer.load_material(
        "test_material", sim, os.path.join(file_dir, "test_data", "materials.yml")
    )

    test_body.material = test_material
    test_body.initial_condition = test_initial_condtion
    test_body.equation = test_eqn

    # detect boundary
    line = get_boundaries_in_box(0, 0, 0, 1, 0, 0, 2, test_domain)
    ph_test_boundary = add_physical_group(1, [line], "test_boundary")

    # initialize test boundary
    test_boundary = elmer.Boundary(sim, "test_boundary", [ph_test_boundary])

    # setup component object
    test_component = elmer.Component(
        sim, "test_component", [test_body], [test_boundary], {"test_parameter": "1.0"}
    )
    test_boundary.data.update({"test_parameter": "1.0"})

    with tempfile.TemporaryDirectory() as tempdir:
        sim.write_startinfo(tempdir)
        sim.write_sif(tempdir)

        files = os.listdir(tempdir)
        # Ensure that we have generated the same files as stored in git.
        assert os.listdir(sim_dir) == files
        for file in files:
            # Check that the contents of each file agrees with the reference.
            assert filecmp.cmp(os.path.join(sim_dir, file), os.path.join(tempdir, file))

    # check format of sif file
    simulation = {}
    constants = {}
    equations = {}
    solvers = {}
    materials = {}
    bodies = {}
    boundaries = {}
    body_forces = {}
    components = {}
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
        "Component",
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
        components,
        initial_conditions,
    ]

    # We already know that the reference sif file is equal to that which we generated,
    # so read the reference sif file here.
    with open(os.path.join(sim_dir, "case.sif"), "r") as f:
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
                        str.strip, re.split(r"(\d+)", line)
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

    assert simulation == sim.settings
    assert equations["test_equation"] == test_eqn.get_data()
    assert solvers["test_solver"] == test_solver.get_data()
    assert solvers["test_post_solver"] == test_post_solver.get_data()
    assert materials["test_material"] == test_material.get_data()
    assert bodies["test_body"] == test_body.get_data()
    assert boundaries["test_boundary"] == test_boundary.get_data()
    assert body_forces["test_body_force"] == test_body_force.get_data()
    assert components["test_component"] == test_component.get_data()
    assert (
        initial_conditions["test_initial_condition"] == test_initial_condtion.get_data()
    )
