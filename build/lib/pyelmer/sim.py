# created by Arved Enders-Seidlitz on 31.07.2020
#
# Run elmer simulations & post processing.

import yaml
import subprocess
import multiprocessing
import platform
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import postprocessing

# settings
simulations_file = './simdata/simulations.yml'


def simulation(simulations_file):
    with open(simulations_file, 'r') as f:
        simulations = yaml.safe_load(f)

    count = multiprocessing.cpu_count()
    if platform.system() == 'Windows':
        count -= 1
    print('Working on ', count, ' cores.')
    pool = multiprocessing.Pool(processes=count)
    pool.map(run_elmer, simulations.values())


def run_elmer(sim_path):
    print('Starting simulation ', sim_path, ' ...')
    with open(sim_path + '/elmersolver.log', 'w') as fp:
        args = ['ElmerSolver.exe', './case.sif']
        subprocess.run(args, cwd=sim_path, stdout=fp, stderr=fp)
    print('Finished simulation ', sim_path, ' .')
    postprocessing.probes(sim_path)
    postprocessing.boundary_scalars(sim_path)
    postprocessing.heat_flux(sim_path)


if __name__ == "__main__":
    simulation(simulations_file)
