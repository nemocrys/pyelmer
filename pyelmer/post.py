import numpy as np
import matplotlib.pyplot as plt
import re
from dataclasses import dataclass

@dataclass
class LinearIteration:
    lines: list


@dataclass
class NonlinearIteration:
    linear_iterations: list
    nrm: float = 0
    relc: float = 0


@dataclass
class SteadyStateIteration:
    nonlinear_iterations: list
    nrm: float = 0
    relc: float = 0

    def nonlin_relc(self):
        relcs = []
        for itr in self.nonlinear_iterations:
            relcs.append(itr.relc)
        return relcs


@dataclass
class SolverResiduals:
    steady_state_iterations: list
    
    def ss_relc(self):
        relcs = []
        for itr in self.steady_state_iterations:
            relcs.append(itr.relc)
        return relcs


def scan_logfile(sim_dir):
    """Scan log file for errors and warnings.

    Args:
        sim_dir (str): Simulation directory

    Returns:
        list[str], list[str]: error messages, warnings
    """
    with open(sim_dir + '/elmersolver.log', 'r') as f:
        log = f.readlines()
    err = []
    warn = []
    for line in log:
        if 'ERROR' in line:
            err.append(line[:-1])
        if 'WARNING' in line:
            err.append(line[:-1])
    return err, warn


def plot_residuals(sim_dir, solvers):
    """Plot residuals in log file

    Args:
        sim_dir (str): Simulation directory
    """
    with open(sim_dir + '/elmersolver.log', 'r') as f:
        log = f.readlines()
    for line in log:
        line = line[:-1]
    
    residuals = {}
    for solver in solvers:
        residuals.update({solver: SolverResiduals([SteadyStateIteration([])])})

    for line in log:
        if 'ComputeChange' in line:
            txt = re.sub(' +', ' ', line.split(' ( ')[-1].split(' ) ')[0]).lstrip().split (' ')
            nrm = float(txt[0])
            relc = float(txt[1])
            for solver in solvers:
                if solver.lower() in line:
                    if 'NS' in line:
                        nli = NonlinearIteration([], nrm, relc)
                        residuals[solver].steady_state_iterations[-1].nonlinear_iterations.append(nli)
                    if 'SS' in line:
                        residuals[solver].steady_state_iterations[-1].nrm = nrm
                        residuals[solver].steady_state_iterations[-1].relc = relc
                        residuals[solver].steady_state_iterations.append(SteadyStateIteration([]))
    
    # remove empty last SteadyStateIteration
    for solver in solvers:
        if residuals[solver].steady_state_iterations[-1].nonlinear_iterations == []:
            del residuals[solver].steady_state_iterations[-1]
        else:
            raise ValueError('Incomplete Steady State iteration in logfile.')

    
    fig, ax = plt.subplots(1, 2)
    for solver in solvers:
        res = residuals[solver]
        # steady state iterations
        line, = ax[0].plot(res.ss_relc(), '-x')
        line.set_label(solver)
        # nonlinear iterations
        for i, ss_iter in enumerate(res.steady_state_iterations):
            line, = ax[1].plot(ss_iter.nonlin_relc(), '-x')
            line.set_label(solver + ' ss ' + str(i))

    ax[0].set_title('Steady State Iterations')
    ax[0].legend()
    ax[0].set_yscale('log')

    ax[1].set_title('Nonlinear Iterations')
    ax[1].legend()
    ax[1].set_yscale('log')
    plt.show()

if __name__ == "__main__":
    sim_dir = r'C:\Users\enders-seidlitz\Documents\GitHub\sim-elmerthermo\simdata\2020-09-21_14-01_vacuum-1st'
    sim_dir = r'C:\Users\enders-seidlitz\Documents\GitHub\sim-elmerthermo\simdata\2020-09-21_16-45_vacuum-2nd_modparams'
    plot_residuals(sim_dir, ['heat equation', 'statmagsolver'])
