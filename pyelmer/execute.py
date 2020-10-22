import os
import shutil
import subprocess
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
        if os.path.exists('C:/Program Files/Elmer 8.4-Release/bin/ElmerGrid.exe'):
            elmergrid = 'C:/Program Files/Elmer 8.4-Release/bin/ElmerGrid.exe'
        else:
            elmergrid = 'ElmerGrid'

    args = [elmergrid, '14', '2', meshfile]
    with open(sim_dir + '/elmergrid.log', 'w') as f:
        subprocess.run(args, cwd=sim_dir, stdout=f, stderr=f)

    mesh_dir = sim_dir + '/' + meshfile.split('.')[0]
    files = os.listdir(mesh_dir)
    for f in files:
        if os.path.exists(sim_dir + '/' + f):
            os.remove(sim_dir + '/' + f)
        shutil.move(mesh_dir + '/' + f, sim_dir)
    shutil.rmtree(mesh_dir)


def run_elmer_solver(sim_dir, elmersolver=None):
    """Run ElmerSolver with input file case.sif

    Args:
        sim_dir (str): Simulation directory
        elmersolver (str): ElmerSolver executable
    """
    if elmersolver is None:
        # On Windows ElmerSolver.exe is not found once gmsh.initialize() was executed.
        # Try to use abs-path instead.
        if os.path.exists('C:/Program Files/Elmer 8.4-Release/bin/ElmerSolver.exe'):
            elmersolver = 'C:/Program Files/Elmer 8.4-Release/bin/ElmerSolver.exe'
        else:
            elmersolver = 'ElmerSolver'

    args = [elmersolver, 'case.sif']
    with open(sim_dir + '/elmersolver.log', 'w') as f:
        subprocess.run(args, cwd=sim_dir, stdout=f, stderr=f)
