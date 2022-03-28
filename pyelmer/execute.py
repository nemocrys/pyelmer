"""Utility functions for the execution of ElmerSolver and ElmerGrid."""
import os
import shutil
import subprocess
import multiprocessing
import platform


def run_elmer_grid(sim_dir, meshfile, elmergrid=None):
    """Run ElmerGrid on gmsh meshfile and move everithing into main
    directory.

    Args:
        sim_dir (str): Simulation directory
        meshfile (str): Filename of .msh file
        elmergrid (str): ElmerGrid executable
    """
    if elmergrid is None:
        # On Windows ElmerGrid.exe is not found once gmsh.initialize() was executed.
        # Try to use abs-path instead.
        if os.path.exists("C:/Program Files/Elmer 9.0-Release/bin/ElmerGrid.exe"):
            elmergrid = "C:/Program Files/Elmer 9.0-Release/bin/ElmerGrid.exe"
        else:
            elmergrid = "ElmerGrid"

    args = [elmergrid, "14", "2", meshfile]
    with open(sim_dir + "/elmergrid.log", "w") as f:
        subprocess.run(args, cwd=sim_dir, stdout=f, stderr=f)

    mesh_dir = sim_dir + "/" + ".".join(meshfile.split(".")[:-1])
    files = os.listdir(mesh_dir)
    for f in files:
        if os.path.exists(sim_dir + "/" + f):
            os.remove(sim_dir + "/" + f)
        shutil.move(mesh_dir + "/" + f, sim_dir)
    shutil.rmtree(mesh_dir)


def run_elmerf90(userfile_in, userfile_out, elmerf90=None):
    """Compile user defined function  or solver in .F90 format to share object (.so) file in linux.
    using elmerf90 compiler
    It has been tested to work in Ubuntu 16.04 computer

    Args:
        userfile_in (str) : Filename of .F90 file
        userfile_out (str) : Filename of .so file
        elmerf90 (str): elmerf90 executable
    """
    if elmerf90 is None:
        # On Windows ElmerSolver.exe is not found once gmsh.initialize() was executed.
        # Try to use abs-path instead.
        if os.path.exists("C:/Program Files/Elmer 9.0-Release/bin/elmerf90.bat"):
            elmerf90 = "C:/Program Files/Elmer 9.0-Release/bin/elmerf90.bat"
        else:
            elmerf90 = "elmerf90"
    subprocess.call([elmerf90, "-o", userfile_out, userfile_in])


def run_elmer_solver(sim_dir, elmersolver=None):
    """Run ElmerSolver with input file case.sif.

    Args:
        sim_dir (str): Simulation directory
        elmersolver (str): ElmerSolver executable
    """
    if elmersolver is None:
        # On Windows ElmerSolver.exe is not found once gmsh.initialize() was executed.
        # Try to use abs-path instead.
        if os.path.exists("C:/Program Files/Elmer 9.0-Release/bin/ElmerSolver.exe"):
            elmersolver = "C:/Program Files/Elmer 9.0-Release/bin/ElmerSolver.exe"
        else:
            elmersolver = "ElmerSolver"

    args = [elmersolver, "case.sif"]
    with open(sim_dir + "/elmersolver.log", "w") as f:
        subprocess.run(args, cwd=sim_dir, stdout=f, stderr=f)


def run_multicore(count, sim_dirs, meshfiles, elmergrid=None, elmersolver=None):
    """Run multiple instances of ElmerGrid, ElmerSolver.

    Args:
        count (int): Number of processes
        sim_dirs (list): Simulation directories
        meshfiles (list): Mesh files for ElmerGrid
        elmergrid (str, optional): Path to executable. Defaults to None.
        elmersolver (str, optional): Path to executable. Defaults to
        None.
    """
    pool = multiprocessing.Pool(processes=count)
    args = []
    for i in range(len(sim_dirs)):
        args.append((sim_dirs[i], meshfiles[i], elmergrid, elmersolver))
    pool.starmap(_run_grid_solver, args)


def _run_grid_solver(sim_dir, meshfile, elmergrid, elmersolver):
    print("Starting ElmerGrid for", meshfile)
    run_elmer_grid(sim_dir, meshfile, elmergrid)
    print("Starting ElmerSolver in", sim_dir)
    run_elmer_solver(sim_dir, elmersolver)
    print("Finished ElmerSolver in", sim_dir)
