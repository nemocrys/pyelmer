import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass


@dataclass
class LinearIteration:
    idx: list
    relc: list


@dataclass
class NonlinearIteration:
    nrm: float = 0
    relc: float = 0
    linear_iteration: LinearIteration = LinearIteration([], [])


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

    def nonlin_nrm(self):
        nrms = []
        for itr in self.nonlinear_iterations:
            nrms.append(itr.nrm)
        return nrms


@dataclass
class SolverResiduals:
    steady_state_iterations: list

    def ss_relc(self):
        relcs = []
        for itr in self.steady_state_iterations:
            relcs.append(itr.relc)
        return relcs

    def ss_nrm(self):
        nrms = []
        for itr in self.steady_state_iterations:
            nrms.append(itr.nrm)
        return nrms


def scan_logfile(sim_dir):
    """Scan log file for errors and warnings.

    Args:
        sim_dir (str): Simulation directory

    Returns:
        list[str], list[str], dict: error messages, warnings, statistics
    """
    with open(sim_dir + "/elmersolver.log", "r") as f:
        log = f.readlines()
    for i in range(len(log)):
        log[i] = log[i][:-1]
    err = []
    warn = []
    stats = {}
    for line in log:
        if "ERROR" in line and line not in err:
            err.append(line)
        if "WARNING" in line and line not in warn:
            warn.append(line)
        if "SOLVER TOTAL TIME(CPU,REAL):" in line:
            s = " ".join(line.split()).split(" ")
            cpu_time = float(s[3])
            real_time = float(s[4])
            stats.update({"CPU-time": cpu_time, "real-time": real_time})
    return err, warn, stats


def plot_residuals(sim_dir, solvers, save=False):
    """Plot residuals in log file. Does not work very well.

    Args:
        sim_dir (str): Simulation directory
        solvers (list): solvers to analyze - currently works only for
        'heat equation' and 'statmagsolver'
    """
    with open(sim_dir + "/elmersolver.log", "r") as f:
        log = f.readlines()
    for i in range(len(log)):
        log[i] = log[i][:-1]

    residuals = {}
    linres_str = {}  # string to detect linear solver residuals
    lin_res_idx = []
    lin_res_relc = []
    for solver in solvers:
        residuals.update({solver: SolverResiduals([SteadyStateIteration([])])})
        if solver == "heat equation":
            linres_str.update({solver: "HeatSolve: Assembly done"})
        if solver == "statmagsolver":
            linres_str.update({solver: "StatMagSolve: Set boundaries done"})

    linres_solver = ""  # currently in section with residuals of this solver
    for line in log:
        if "ComputeChange" in line:
            txt = (
                re.sub(" +", " ", line.split(" ( ")[-1].split(" ) ")[0])
                .lstrip()
                .split(" ")
            )
            nrm = float(txt[0])
            relc = float(txt[1])
            for solver in solvers:
                if solver.lower() in line:
                    if "NS" in line:
                        nli = NonlinearIteration(nrm, relc)
                        residuals[solver].steady_state_iterations[
                            -1
                        ].nonlinear_iterations.append(nli)
                    if "SS" in line:
                        residuals[solver].steady_state_iterations[-1].nrm = nrm
                        residuals[solver].steady_state_iterations[-1].relc = relc
                        residuals[solver].steady_state_iterations.append(
                            SteadyStateIteration([])
                        )
            if linres_solver:
                li = LinearIteration(lin_res_idx, lin_res_relc)
                residuals[linres_solver].steady_state_iterations[
                    -1
                ].nonlinear_iterations[-1].linear_iteration = li
                linres_solver = ""
        else:
            for solver in solvers:
                if line == linres_str[solver]:
                    linres_solver = solver
                    lin_res_idx = []
                    lin_res_relc = []
                elif linres_solver == solver:
                    txt = line.lstrip().split(" ")
                    try:
                        lin_res_idx.append(float(txt[0]))
                        lin_res_relc.append(float(txt[1]))
                    except ValueError:
                        print("Problem in evaluation of linear residuals.")
                        print("Could not read the following line:\n", line)

    # remove empty last SteadyStateIteration
    for solver in solvers:
        if residuals[solver].steady_state_iterations[-1].nonlinear_iterations == []:
            del residuals[solver].steady_state_iterations[-1]
        else:
            print("WARNING: Incomplete Steady State iteration in logfile.")

    figs = []
    axes = []

    fig_ss, ax_ss = plt.subplots(1, 2, figsize=(5.75, 4))
    figs.append(fig_ss)
    axes.append(ax_ss)
    # ax_ss.set_title('Steady State Iterations')
    fig_ss.title = "residuals_ss"

    fig_ns, ax_ns = plt.subplots(1, 2, figsize=(5.75, 4))
    figs.append(fig_ns)
    axes.append(ax_ns)
    # ax_ns.set_title('Nonlinear Iterations')
    fig_ns.title = "residuals_ns"

    fig_li, ax_li = plt.subplots(1, len(linres_str), figsize=(5.75, 4))
    for i, solver in enumerate(linres_str):
        ax_li[i].set_title(solver)
    figs.append(fig_li)
    fig_li.title = "residuals_li"

    for solver in solvers:
        res = residuals[solver]
        # steady state iterations
        (line,) = ax_ss[0].plot(res.ss_relc(), "-x")
        line.set_label(solver)
        (line,) = ax_ss[1].plot(res.ss_nrm(), "-x")
        line.set_label(solver)
        # nonlinear iterations
        for i, ss_iter in enumerate(res.steady_state_iterations):
            (line,) = ax_ns[0].plot(ss_iter.nonlin_relc(), "-x")
            line.set_label(solver + " ss " + str(i))
            (line,) = ax_ns[1].plot(ss_iter.nonlin_nrm(), "-x")
            line.set_label(solver + " ss " + str(i))

    for i, solver in enumerate(linres_str):
        for i_ss, ss_iter in enumerate(residuals[solver].steady_state_iterations):
            for i_nl, nl_iter in enumerate(ss_iter.nonlinear_iterations):
                li = nl_iter.linear_iteration
                (line,) = ax_li[i].plot(li.idx, li.relc)
                line.set_label("ss " + str(i_ss) + " nl " + str(i_nl))

    for ax in axes:
        ax[1].legend()
        ax[0].set_yscale("log")
        ax[0].set_xlabel("iteration")
        ax[1].set_xlabel("iteration")
        ax[0].set_ylabel("RELC")
        ax[1].set_ylabel("NRM")
    for ax in ax_li:
        # ax.legend()
        ax.set_yscale("log")
        ax.set_xlabel("iteration")
        ax.set_ylabel("RELC")

    fig_dir = sim_dir + "/results/"
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    for fig in figs:
        fig.tight_layout()
        if save:
            fig.savefig(fig_dir + fig.title + ".pdf")
            fig.savefig(fig_dir + fig.title + ".png")

    return figs, axes


def dat_to_dataframe(dat_file):
    """Read a .dat file generated by elmer (e.g. SaveData/SaveScalars
    solver) into a pandas dataframe.

    Args:
        dat_file (str): file path of elmer .dat file

    Returns:
        pandas.core.frame.DataFrame: content of dat file with header
    """
    with open(f"{dat_file}.names") as f:
        lines = f.readlines()
    names = []
    names_start = False
    for line in lines:
        if names_start == True:
            names.append(line.split(":")[-1].strip())
        if "Data on different columns" in line:
            names_start = True
    return pd.read_table(dat_file, names=names, delim_whitespace=True)
