# created by Arved Enders-Seidlitz on 14.07.2020
#
# Mapping of Elmer simulation input file to python.


class Simulation:
    def __init__(self):
        self.materials = {}
        self.bodies = {}
        self.boundaries = {}
        self.body_forces = {}
        self.initial_conditions = {}
        self.solvers = {}
        self.equations = {}
        self.constants = {
         'Stefan Boltzmann': 5.6704e-08
        }
        self.simulation_settings = {
            "Max Output Level": 4,
            "Coordinate System": "Axi Symmetric",
            "Simulation Type": "Steady state",
            "Steady State Max Iterations": 3,
            "Timestepping Method": "BDF",
            "BDF Order" : 2,
            "Output Intervals" : 1,
            "Timestep Intervals": 10,
            "Timestep Sizes": 100
        }

    def write_sif(self, simulation_dir, filename):
        self._set_ids()
        with open(simulation_dir + '/' + filename, 'w') as f:
            f.write('''Header\n  CHECK KEYWORDS "Warn"\n  Mesh DB "." "."\nEnd\n\n''')
            f.write('Simulation\n')
            f.write(self._dict_to_str(self.simulation_settings))
            f.write('End\n\n')
            f.write('Constants\n')
            f.write(self._dict_to_str(self.constants))
            f.write('End\n\n')
            for equation_name, equation in self.equations.items():
                f.write('! ' + equation_name + '\n')
                f.write('Equation ' + str(equation.id) + '\n')
                f.write(self._dict_to_str(equation.get_data()))
                f.write('End\n\n')
            f.write('\n')
            for solver_name, solver in self.solvers.items():
                f.write('! ' + solver_name + '\n')
                f.write('Solver ' + str(solver.id) + '\n')
                f.write(self._dict_to_str(solver.data))
                f.write('End\n\n')
            f.write('\n')
            for material_name, material in self.materials.items():
                f.write('! ' + material_name + '\n')
                f.write('Material ' + str(material.id) + '\n')
                f.write(self._dict_to_str(material.data))
                f.write('End\n\n')
            f.write('\n')
            for body_name, body in self.bodies.items():
                f.write('! ' + body_name + '\n')
                f.write('Body ' + str(body.id) + '\n')
                f.write(self._dict_to_str(body.get_data()))
                f.write('End\n\n')
            f.write('\n')
            for boundary_name, boundary in self.boundaries.items():
                f.write('! ' + boundary_name + '\n')
                f.write('Boundary Condition ' + str(boundary.id) + '\n')
                f.write(self._dict_to_str(boundary.get_data()))
                f.write('End\n\n')
            f.write('\n')
            for body_force_name, body_force in self.body_forces.items():
                f.write('! ' + body_force_name + '\n')
                f.write('Body Force ' + str(body_force.id) + '\n')
                f.write(self._dict_to_str(body_force.get_data()))
                f.write('End\n\n')
            f.write('\n')
            for initial_condition_name, initial_condition in self.initial_conditions.items():
                f.write('! ' + initial_condition_name + '\n')
                f.write('Initial Condition ' + str(initial_condition.id) + '\n')
                f.write(self._dict_to_str(initial_condition.data))
                f.write('End\n\n')
        print('Wrote ', filename)

    def write_boundaries(self, simulation_dir, filename):
        with open(simulation_dir + '/' + filename, 'w') as f:
            f.write('ID    Boundary-Name\n')
            for key, boundary in self.boundaries.items():
                f.write(''.join([str(boundary.id), '    ', key, '\n']))

    def _dict_to_str(self, dictionary):
        text = ''
        for key, value in dictionary.items():
            text = ''.join([text, '  ', key, ' = ', str(value), '\n'])
        return text
    
    def _set_ids(self):
        objects = [self.solvers, self.equations, self.materials, self.bodies, self.boundaries,
                   self.body_forces, self.initial_conditions]
        # give each object id
        for obj in objects:
            id = 1
            for key in obj:
                obj[key].id = id
                id += 1


class Body:
    def __init__(self, simulation, name, body_id=0):
        self.id = 0
        self.body_id = body_id
        self.equation = None
        self.initial_condition = None
        self.material = None
        self.body_force = None
        self.name = name
        simulation.bodies.update({name: self})

    def get_data(self):
        d = {'Target Bodies(1)': self.body_id}
        if not self.equation is None:
            d.update({'Equation': self.equation.id})
        if not self.initial_condition  is None:
            d.update({'Initial Condition': self.initial_condition.id})
        if not self.material is None:
            d.update({'Material': self.material.id})
        if not self.body_force is None:
            d.update({'Body Force': self.body_force.id})
        return d


class Boundary:
    def __init__(self, simulation, name, surf_id=0):
        self.id = 0
        self.surface_id = surf_id
        #self.material = None
        self.radiation = False
        self.fixed_temperature = None
        self.zero_potential = False
        self.save_flux = True
        self.smart_heater = False
        self.smart_heater_T = 0
        self.phase_change = False
        self.phase_change_vel = 0
        self.material = None
        self.normal_target_body = None
        self.phase_change_body = None
        self.name = name
        simulation.boundaries.update({name: self})
    
    def get_data(self):
        d = {}
        d.update({'Target Boundaries(1)': self.surface_id})
        if self.radiation:
            d.update({'Heat Flux BC': 'True'})
            d.update({'Radiation': 'Diffuse Gray'})
            # use definition of emissivity in material
            # d.update({'Radiation Target Body': -1})
            # d.update({'Emissivity': self.material.data['Emissivity']})
        if not self.fixed_temperature is None:
            d.update({'Temperature': self.fixed_temperature})
        if self.zero_potential:
            d.update({'Potential 1': 0, 'Potential 2': 0})
        if self.save_flux:
            d.update({'Save Scalars': 'Logical True', 'Save Line': 'Logical True'})
        if self.smart_heater:
            d.update({'Smart Heater Boundary': 'Logical True',
                      'Smart Heater Temperature': self.smart_heater_T})
        if self.phase_change:
            d.update({
                'Phase Change': 'Logical True',
                'Phase Velocity 1': 0,
                'Phase Velocity 2': self.phase_change_vel,
                'Melting Point': self.material.data['Melting Point'],
                'Latent Heat': self.material.data['Latent Heat'],
                'Normal Target Body': self.normal_target_body.id,
                'Heat Flux': 'Variable Coordinate 1\n    Real Procedure "SteadyPhaseChange" "MeltingHeat"',
                'Body Id': 'Integer ' + str(self.phase_change_body.id)
            })
        return d


class Material:
    def __init__(self, simulation, name, data={}):
        self.id = 0
        self.data = data
        self.name = name
        simulation.materials.update({name: self})


class BodyForce:
    def __init__(self, simulation, name, data={}):
        self.id = 0
        self.joule_heat = False
        self.current_density = 0
        self.heat_source = 0
        self.integral_heat_source = 0
        self.smart_heat_control = False
        self.smart_heater_control_point = []
        self.smart_heater_T = 0
        self.name = name
        simulation.body_forces.update({name: self})
    def get_data(self):
        d = {}
        if self.joule_heat:
            d.update({'Joule Heat': 'Logical True'})
        if self.current_density != 0:
            d.update({'Current Density': self.current_density})
        if self.heat_source != 0:
            d.update({'Heat Source': self.heat_source})
        if self.integral_heat_source != 0:
            d.update({'Integral Heat Source': self.integral_heat_source})
        if self.smart_heat_control:
            d.update({'Smart Heater Control': 'Logical True'})
            if self.joule_heat:
                d.update({'Heat Source': '0'})
        if self.smart_heater_control_point != []:
            cp = self.smart_heater_control_point
            d.update({'Smart Heater Control Point(3)': str(cp[0]) + ' ' + str(cp[1]) + ' ' + str(cp[2]),
                      'Smart Heater Temperature': self.smart_heater_T})
        return d


class InitialCondition:
    def __init__(self, simulation, name, data={}):
        self.id = 0
        self.data = data
        self.name = name
        simulation.initial_conditions.update({name: self})


class Solver:
    def __init__(self, simulation, name, data={}):
        self.id = 0
        self.data = data
        self.name = name
        simulation.solvers.update({name: self})


class Equation:
    def __init__(self, simulation, name, solvers):
        self.id = 0
        self.solvers = solvers
        self.name = name
        simulation.equations.update({name: self})

    def get_data(self):
        n_solvers = len(self.solvers)
        solver_id_str = ''
        for solver in self.solvers:
            solver_id_str += str(solver.id) + ' '
        return {'Active Solvers(' + str(n_solvers) + ')': solver_id_str}
