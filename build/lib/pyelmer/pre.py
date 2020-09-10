# created by Arved Enders-Seidlitz on 31.07.2020
#
# Create elmer simulations.

import os
import shutil
import datetime
import subprocess
import copy
import yaml
import ruamel.yaml  # because pyYAML doesn't support numbers of style 1e3
# import cz_2d


###############################
# settings
simulation_file = './simulation.yml'
base_params_file = './base_parameters.yml'
###############################


def create_simulations(simulation_file, base_params_file):
    """Generate input for Elmer simulation of 2D Czochalski growth for case specified
    in simulation_file, using the parameters in base_params_file.
    If simulation_file containts lists, permutations of the listed values are generated
    and multiple simulations 'parameter studies' are created.

    Args:
        simulation_file (str): file path to yaml file
        base_params_file (str): file path to yaml file
    """
    with open(simulation_file, 'r') as f:
        sim_dict = ruamel.yaml.safe_load(f)
    with open(base_params_file, 'r') as f:
        base_dict = ruamel.yaml.safe_load(f)

    if not os.path.exists('./simdata'):
        os.mkdir('./simdata')

    # check if there are lists -> 'parameter study' case
    param_lists_dict = search_parameter_lists(sim_dict)

    if param_lists_dict != {}:  # parameter study
        simulations = create_permutations(param_lists_dict, base_dict, sim_dict)
        base_path = './simdata/' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M") + '_ps_' + sim_dict['settings']['name'] +'/'
        os.mkdir(base_path)
        shutil.copy2('./simulation.yml', base_path + 'simulation.yml')
    else:  # single simulation
        simulation = replace_params(base_dict, sim_dict)
        simulations = [simulation]
        base_path = './simdata/' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")  + '_'

    sim_dict = {}
    for sim in simulations:
        sim_path = base_path + sim['settings']['name'] + '/'
        os.mkdir(sim_path)

        sim['settings']['path'] = sim_path
        parameters_file = sim_path + 'parameters.yml'
        with open(parameters_file, 'w') as f:
            yaml.dump(sim, f, sort_keys=False)

        # begin workaround:
        # execute cz_2d.py using subprocess to avoid problems caused by gmsh.initialize()
        # (ElmerGrid.exe is not found anymore)
        args = ['python', './cz_2d.py', parameters_file]
        subprocess.run(args)
        mesh_file = sim_path + sim['settings']['name'] + '.msh'
        # once the bug is fixed use instead:
        # mesh_file = cz_2d.geometry_sif(parameters_file)
        # end workaround

        with open(sim_path + 'elmergrid.log', 'w') as f:
            args = ['ElmerGrid.exe', '14', '2', mesh_file]
            subprocess.run(args, stdout=f, stderr=f)
        # move everything in main directory
        mesh_folder = sim_path + sim['settings']['name'] + '/'
        files = os.listdir(mesh_folder)
        for f in files:
            shutil.move(mesh_folder+f, sim_path)
        shutil.rmtree(mesh_folder)

        shutil.copy2('./cz_2d.py', sim_path + 'cz_2d.py')
        shutil.copy2('./ELMERSOLVER_STARTINFO', sim_path + 'ELMERSOLVER_STARTINFO')

        sim_dict.update({sim['settings']['name']: sim_path})

    # export simulations.yml
    with open('./simdata/simulations.yml', 'w') as f:
        yaml.dump(sim_dict, f, sort_keys=False)


def replace_params(base_dict, sim_dict):
    """Replaces entries of base_dict that are listed in sim_dict with
    the values given in sim_dict. Recursive function to cover the whole
    depth of the dict.

    Args:
        base_dict (dict): Dictionary containing all simulation parameters
        sim_dict (dict): Dictionary containing the parameters differing from base

    Raises:
        ValueError: If a key of sim_dict is not contained in base_dict

    Returns:
        dict: Modified version of base_dict, now containing the values of sim_dict
        at respective keys.
    """
    for key in sim_dict:
        if not key in base_dict.keys():
            raise ValueError('Wrong parameter name in simulation.yml: ' + key)
        if type(sim_dict[key]) is dict:
            base_dict[key] = replace_params(base_dict[key], sim_dict[key])
        else:
            base_dict[key] = sim_dict[key]
    return base_dict


def search_parameter_lists(sim_dict):
    """Searches for lists in values of sim dict. Only key-value pairs
    with lists are kept. Recursive function to cover the whole depth of
    the dict.

    Args:
        sim_dict (dict): Dictionary containing the simulation parameters

    Returns:
        dict: dict (of dicts) with only list as values.
    """
    param_lists_dict = {}
    for key in sim_dict:
        if type(sim_dict[key]) is dict:
            params = search_parameter_lists(sim_dict[key])
            if params != {}:
                param_lists_dict.update({key: params})
        if type(sim_dict[key]) is list:
            param_lists_dict.update({key: sim_dict[key]})
    return param_lists_dict


def get_first_parameter(sim_dict):
    """For the setup of the base simulation for the permutations,
    a set of parameters containing the first value of each list
    is required. For the whole depth of the dict, the following
    is done: If value = dict - recursion, if value = list - take
    first element, else - keep entry as is.

    Args:
        sim_dict (dict): Dictionary containing simulation parameters

    Returns:
        dict: Simulation parameter dict without lists
    """
    base_params = {}
    for key in sim_dict:
        if type(sim_dict[key]) is dict:
            params = get_first_parameter(sim_dict[key])
            if params != {}:
                base_params.update({key: params})
        elif not type(sim_dict[key]) is list:
            base_params.update({key: sim_dict[key]})
    return base_params


def create_variation_lists(param_lists_dict):
    # TODO dirty solution! -> recursive?
    variations = []
    for key_0, value_0 in param_lists_dict.items():
        if type(value_0) is list:
            variation = []
            for element in value_0:
                variation.append({key_0: element})
            variations.append(variation)
        else:
            for key_1, value_1 in value_0.items():
                if type(value_1) is list:
                    variation = []
                    for element in value_1:
                        variation.append({key_0: {key_1: element}})
                    variations.append(variation)
                else:
                    for key_2, value_2 in value_1.items():
                        if type(value_2) is list:
                            variation = []
                            for element in value_2:
                                variation.append({key_0: {key_1: {key_2: element}}})
                            variations.append(variation)
                        else:
                            raise ValueError('To high depth in data dict.')
    return variations


def create_permutations(param_lists_dict, base_dict, sim_dict):
    """Creates permutations of all listed parameters.

    Args:
        param_lists_dict (dict): Dict containing parameter lists for respective keys.
        base_dict (dict): Dict containing base parameters.
        sim_dict (dict): Dict containing user defined parameters for this simulation (both lists and values only).

    Returns:
        list: Dict with all parameter permutations.
    """
    single_params = get_first_parameter(sim_dict)
    base = replace_params(base_dict, single_params)

    variations = create_variation_lists(param_lists_dict)
    
    simulations = [base]
    for variation in variations:
        simulations_new = []
        for sim in simulations:
            for variant in variation:
                sim_new = replace_params(copy.deepcopy(sim), variant)
                print(variant)
                name = str(variant).replace('{','').replace('}','').replace(' ', '').replace("'", '').replace('.', '-').split(':')
                name = name[-2] + '_' + name[-1]
                print(name)
                sim_new['settings']['name'] += '_' + name
                simulations_new.append(sim_new)
        simulations = simulations_new       
    
    return simulations


if __name__ == "__main__":
    create_simulations(simulation_file, base_params_file)
