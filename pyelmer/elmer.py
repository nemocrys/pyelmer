"""Wrapper for Elmer simulation input file (sif-file).

The sif file is represented by a Simulation object, which contains
dictionaries of objects for each section (i.e. Solvers, Bodies, 
Materials, etc.).
Each of these objects contains the data required in the sif in form
of a dictionary called data.
For some objects (e.g. Bodies) this dictionary is written
automatically, provided that the required member variables of the
object are set.
"""

import os
import yaml


DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


class Simulation:
    """Main wrapper for the sif file. The simulation class is used to
    collect all information.
    In the end, it writes the sif-file.
    """

    def __init__(self):
        self.materials = {}
        self.bodies = {}
        self.boundaries = {}
        self.body_forces = {}
        self.initial_conditions = {}
        self.solvers = {}
        self.equations = {}
        self.constants = {"Stefan Boltzmann": 5.6704e-08}
        self.settings = {}

    def write_sif(self, simulation_dir):
        """Write sif file.

        Args:
            simulation_dir (str): Path of simulation directory
        """
        self._set_ids()
        with open(simulation_dir + "/case.sif", "w") as f:
            f.write("""Header\n  CHECK KEYWORDS "Warn"\n  Mesh DB "." "."\nEnd\n\n""")
            f.write("Simulation\n")
            f.write(self._dict_to_str(self.settings))
            f.write("End\n\n")
            f.write("Constants\n")
            f.write(self._dict_to_str(self.constants))
            f.write("End\n\n")
            for equation_name, equation in self.equations.items():
                f.write("! " + equation_name + "\n")
                f.write("Equation " + str(equation.id) + "\n")
                f.write(self._dict_to_str(equation.get_data()))
                f.write("End\n\n")
            f.write("\n")
            for solver_name, solver in self.solvers.items():
                f.write("! " + solver_name + "\n")
                f.write("Solver " + str(solver.id) + "\n")
                f.write(self._dict_to_str(solver.get_data()))
                f.write("End\n\n")
            f.write("\n")
            for material_name, material in self.materials.items():
                f.write("! " + material_name + "\n")
                f.write("Material " + str(material.id) + "\n")
                f.write(self._dict_to_str(material.get_data()))
                f.write("End\n\n")
            f.write("\n")
            for body_name, body in self.bodies.items():
                f.write("! " + body_name + "\n")
                f.write("Body " + str(body.id) + "\n")
                f.write(self._dict_to_str(body.get_data()))
                f.write("End\n\n")
            f.write("\n")
            for boundary_name, boundary in self.boundaries.items():
                f.write("! " + boundary_name + "\n")
                f.write("Boundary Condition " + str(boundary.id) + "\n")
                f.write(self._dict_to_str(boundary.get_data()))
                f.write("End\n\n")
            f.write("\n")
            for body_force_name, body_force in self.body_forces.items():
                f.write("! " + body_force_name + "\n")
                f.write("Body Force " + str(body_force.id) + "\n")
                f.write(self._dict_to_str(body_force.get_data()))
                f.write("End\n\n")
            f.write("\n")
            for (
                initial_condition_name,
                initial_condition,
            ) in self.initial_conditions.items():
                f.write("! " + initial_condition_name + "\n")
                f.write("Initial Condition " + str(initial_condition.id) + "\n")
                f.write(self._dict_to_str(initial_condition.get_data()))
                f.write("End\n\n")
        print("Wrote sif-file.")

    def write_startinfo(self, simulation_dir):
        """Write ELMERSOLVER_STARTINFO file in simulation directory.

        Args:
            simulation_dir (str): simulation directory
        """
        with open(simulation_dir + "/ELMERSOLVER_STARTINFO", "w") as f:
            f.write("case.sif\n")

    def write_boundary_ids(self, simulation_dir):
        """Write yaml-file containing the boundary names and the
        assigned elmer body-ids.

        Args:
            simulation_dir (str): Path of simulation directory
        """
        data = {boundary.id: name for name, boundary in self.boundaries.items()}
        with open(simulation_dir + "/boundaries.yml", "w") as f:
            yaml.dump(data, f, sort_keys=False)

    def _dict_to_str(self, dictionary):
        text = ""
        for key, value in dictionary.items():
            text = "".join([text, "  ", key, " = ", str(value), "\n"])
        return text

    def _set_ids(self):
        objects = [
            self.solvers,
            self.equations,
            self.materials,
            self.bodies,
            self.boundaries,
            self.body_forces,
            self.initial_conditions,
        ]
        # give each object id
        for obj in objects:
            id = 1
            for key in obj:
                obj[key].id = id
                id += 1


class Body:
    """Wrapper for bodies in sif-file."""

    def __init__(self, simulation, name, body_ids=[]):
        """Create body object.

        Args:
            simulation (Simulation Object): The body is added to this
                                            simulation object.
            name (str): Name of the body
            body_ids (list of int): Ids of bodies in mesh.
        """
        simulation.bodies.update({name: self})
        self.id = 0
        self.name = name
        self.body_ids = body_ids
        self.data = {}
        # optional parameters
        self.equation = None
        self.initial_condition = None
        self.material = None
        self.body_force = None

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        d = {
            f"Target Bodies({len(self.body_ids)})": " ".join(
                [str(x) for x in self.body_ids]
            )
        }
        if self.equation is not None:
            d.update({"Equation": f"{self.equation.id}  ! {self.equation.name}"})
        if self.initial_condition is not None:
            d.update(
                {
                    "Initial Condition": f"{self.initial_condition.id}  ! {self.initial_condition.name}"
                }
            )
        if self.material is not None:
            d.update({"Material": f"{self.material.id}  ! {self.material.name}"})
        if self.body_force is not None:
            d.update({"Body Force": f"{self.body_force.id}  ! {self.body_force.name}"})
        d.update(self.data)
        return d


class Boundary:
    """Wrapper for boundaries in sif-file."""

    def __init__(self, simulation, name, geo_ids=[]):
        """Create boundary object.

        Args:
            simulation (Simulation Object): The boundary is added to
                                            this simulation object.
            name (str): Name of the body
            surf_ids (list of int): Ids of boundaries in mesh.
        """
        simulation.boundaries.update({name: self})
        self.id = 0
        self.name = name
        self.data = {}
        # optional parameters
        self.geo_ids = geo_ids
        self.radiation = False
        self.radiation_idealized = False
        self.fixed_temperature = None
        self.fixed_heatflux = None
        self.zero_potential = False
        self.save_scalars = False
        self.save_line = False
        self.smart_heater = False
        self.smart_heater_T = 0
        self.phase_change_steady = False
        self.phase_change_transient = False
        self.phase_change_vel = 0
        self.material = None
        self.normal_target_body = None
        self.phase_change_body = None
        self.heat_transfer_coefficient = 0
        self.T_ext = 0
        self.mesh_update = []

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        d = {
            f"Target Boundaries({len(self.geo_ids)})": " ".join(
                [str(x) for x in self.geo_ids]
            )
        }
        if self.radiation:
            d.update({"Radiation": "Diffuse Gray"})
        if self.radiation_idealized:
            d.update({"Radiation": "Idealized", "External Temperature": self.T_ext})
        if self.fixed_heatflux is not None:
            d.update({"Heat Flux BC": True, "Heat Flux": self.fixed_heatflux})
        if self.fixed_temperature is not None:
            d.update({"Temperature": self.fixed_temperature})
        if self.zero_potential:
            # d.update({'Potential 1': 0, 'Potential 2': 0})
            d.update({"Potential Re": 0, "Potential Im": 0})
        if self.save_scalars:
            d.update({"Save Scalars": "Logical True"})
        if self.save_line:
            d.update({"Save Line": "Logical True"})
        if self.smart_heater:
            d.update(
                {
                    "Smart Heater Boundary": "Logical True",
                    "Smart Heater Temperature": self.smart_heater_T,
                }
            )
        if self.phase_change_steady:
            d.update(
                {
                    "Phase Change": "Logical True",
                    "Phase Velocity 1": 0,
                    "Phase Velocity 2": self.phase_change_vel,
                    "Melting Point": self.material.data["Melting Point"],
                    "Latent Heat": self.material.data["Latent Heat"],
                    "Normal Target Body": self.normal_target_body.id,
                    "Heat Flux": 'Variable Coordinate 1\n    Real Procedure "SteadyPhaseChange" "MeltingHeat"',
                    "Mesh Update 1": 0,
                    "Mesh Update 2": "Equals PhaseSurface",
                    "Body Id": "Integer " + str(self.phase_change_body.id),
                }
            )
        if self.phase_change_transient:
            d.update(
                {
                    # 'Phase Velocity 1': 0,
                    # 'Phase Velocity 2': self.phase_change_vel,
                    "Temperature": self.material.data["Melting Point"],
                    "Normal Target Body": self.normal_target_body.id,
                    # 'Latent Heat': self.material.data['Latent Heat'],
                    # 'Heat Flux': 'Variable Coordinate 1\n    Real Procedure "SteadyPhaseChange" "MeltingHeat"',
                    "Mesh Update 1": 0,
                    "Mesh Update 2": "Equals PhaseSurface",
                    "Body Id": "Integer " + str(self.phase_change_body.id),
                }
            )
        if self.heat_transfer_coefficient != 0:
            d.update(
                {
                    "Heat Transfer Coefficient": self.heat_transfer_coefficient,
                    "External Temperature": self.T_ext,
                }
            )
        if self.mesh_update:
            if len(self.mesh_update) >= 2:
                if self.mesh_update[0] is not None:
                    d.update({"Mesh Update 1": self.mesh_update[0]})
                if self.mesh_update[1] is not None:
                    d.update({"Mesh Update 2": self.mesh_update[1]})
            if len(self.mesh_update) == 3:
                if self.mesh_update[2] is not None:
                    d.update({"Mesh Update 3": self.mesh_update[2]})
        d.update(self.data)
        return d


class Material:
    """Wrapper for materials in sif-file."""

    def __init__(self, simulation, name, data={}):
        """Create material object

        Args:
            simulation (Simulation Object): The material is added to
                                            this simulation object.
            name (str): Name of the material
            data (dict): Material data as in sif-file.
        """
        if simulation is not None:
            simulation.materials.update({name: self})
        self.id = 0
        self.name = name
        self.data = data

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        return self.data


class BodyForce:
    """Wrapper for body forces in sif-file."""

    def __init__(self, simulation, name, data={}):
        """Create body force object.

        Args:
            simulation (Simulation Object): The body force is added to
                                            this simulation object.
            name (str): Name of the body force
            data (dict): Body force data as in sif-file.
        """
        simulation.body_forces.update({name: self})
        self.id = 0
        self.name = name
        self.data = data
        # optional parameters
        self.joule_heat = False
        self.current_density = 0
        self.heat_source = 0
        self.integral_heat_source = 0
        self.smart_heat_control = False
        self.smart_heater_control_point = []
        self.smart_heater_T = 0

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        d = {}
        if self.joule_heat:
            d.update({"Joule Heat": "Logical True"})
        if self.current_density != 0:
            d.update({"Current Density": self.current_density})
        if self.heat_source != 0:
            d.update({"Heat Source": self.heat_source})
        if self.integral_heat_source != 0:
            d.update({"Integral Heat Source": self.integral_heat_source})
        if self.smart_heat_control:
            d.update({"Smart Heater Control": "Logical True"})
            if self.joule_heat:
                d.update({"Heat Source": "0"})
        if self.smart_heater_control_point != []:
            cp = self.smart_heater_control_point
            d.update(
                {
                    "Smart Heater Control Point(3)": str(cp[0])
                    + " "
                    + str(cp[1])
                    + " "
                    + str(cp[2]),
                    "Smart Heater Temperature": self.smart_heater_T,
                }
            )
        d.update(self.data)
        return d


class InitialCondition:
    """Wrapper for initial condition in sif-file."""

    def __init__(self, simulation, name, data={}):
        """Create initial condition object.

        Args:
            simulation (Simulation Object): The initial condition is
                                            added to this simulation
                                            object.
            name (str): Name of the initial condition
            data (dict): Initial condition data as in sif-file.
        """
        simulation.initial_conditions.update({name: self})
        self.id = 0
        self.name = name
        self.data = data

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        return self.data


class Solver:
    """Wrapper for solver in sif-file."""

    def __init__(self, simulation, name, data={}):
        """Create solver object

        Args:
            simulation (Simulation Object): The solver is added to
                                            this simulation object.
            name (str): Name of the solver
            data (dict): Solver data as in sif-file.
        """
        if simulation is not None:
            simulation.solvers.update({name: self})
        self.id = 0
        self.name = name
        self.data = data

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        return self.data


class Equation:
    """Wrapper for equations in sif-file."""

    def __init__(self, simulation, name, solvers):
        """Create equation object

        Args:
            simulation (Simulation Object): The equation is added to
                                            this simulation object.
            name (str): Name of the equation
            solvers ([list]): Solvers in this equation
        """
        simulation.equations.update({name: self})
        self.id = 0
        self.name = name
        self.solvers = solvers
        self.data = {}

    def get_data(self):
        """Generate dictionary with data for sif-file."""
        solver_id_str = ""
        solver_name_str = ""
        for solver in self.solvers:
            solver_id_str += f"{solver.id} "
            solver_name_str += f"{solver.name}, "
        return {
            f"Active Solvers({len(self.solvers)})": f"{solver_id_str}  ! {solver_name_str}"
        }.update(self.data)


def load_simulation(name, setup_file=None):
    """Load simulation settings from database.

    Args:
        name (str): Name of the simulation in database.
        setup_file (str, optional): Path to yaml file cotaining setup.

    Returns:
        Simulation object.
    """
    if setup_file is None:
        setup_file = DATA_DIR + "/simulations.yml"
    with open(setup_file) as f:
        settings = yaml.safe_load(f)[name]
    sim = Simulation()
    sim.settings = settings
    return sim


def load_material(name, simulation=None, setup_file=None):
    """Load material from data base and add it to simulation.

    Args:
        name (str): material name
        simulation (Simulation object)
        setup_file (str, optional): Path to yaml file cotaining setup.

    Returns:
        Material object.
    """
    if setup_file is None:
        setup_file = DATA_DIR + "/materials.yml"
    with open(setup_file) as f:
        data = yaml.safe_load(f)[name]
    return Material(simulation, name, data)


def load_solver(name, simulation=None, setup_file=None):
    """Load solver from data base and add it to simulation.

    Args:
        name (str): solver name
        simulation (Simulation object)
        setup_file (str, optional): Path to yaml file cotaining setup.

    Returns:
        Solver object.
    """
    if setup_file is None:
        setup_file = DATA_DIR + "/solvers.yml"
    with open(setup_file) as f:
        data = yaml.safe_load(f)
        data = data[name]
    return Solver(simulation, name, data)
