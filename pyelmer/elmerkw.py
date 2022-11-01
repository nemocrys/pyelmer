"""Implementation of solver-specific keywords for backward compatibility.
If you don't know what this is please use elmer.py instead."""
from pyelmer import elmer
from pyelmer.elmer import (
    Simulation,
    Solver,
    Body,
    Material,
    Equation,
    InitialCondition,
    StringFromList,
    load_material,
    load_simulation,
    load_solver,
)


DeprecationWarning(
    "Solver-specific keywords will not be developed any further and are just kept for backward compability."
)


class Boundary(elmer.Boundary):
    """Wrapper for boundaries in sif-file."""

    def __init__(self, simulation, name, geo_ids=None):
        """Create boundary object.

        Args:
            simulation (Simulation Object): The boundary is added to
                                            this simulation object.
            name (str): Name of the body
            surf_ids (list of int): Ids of boundaries in mesh.
        """
        super().__init__(simulation, name, geo_ids)
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
        d = super().get_data()
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
        return d


class BodyForce(elmer.BodyForce):
    """Wrapper for body forces in sif-file."""

    def __init__(self, simulation, name, data=None):
        """Create body force object.

        Args:
            simulation (Simulation Object): The body force is added to
                                            this simulation object.
            name (str): Name of the body force
            data (dict): Body force data as in sif-file.
        """
        super().__init__(simulation, name, data)
        self.joule_heat = False
        self.current_density = 0
        self.heat_source = 0
        self.integral_heat_source = 0
        self.smart_heat_control = False
        self.smart_heater_control_point = []
        self.smart_heater_T = 0

    def get_data(self):
        d = super().get_data()
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
        return d
