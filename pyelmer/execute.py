import os
import shutil
import subprocess


def run_elmer_grid(sim_dir, meshfile, elmergrid='win'):
    """Run ElmerGrid on gmsh meshfile and move everithing into main
    directory.

    Args:
        sim_dir (str): Simulation directory
        meshfile (str): Filename of .msh file
        elmergrid (str, optional): ElmerGrid executable
    """
    # On Windows ElmerGrid.exe is not found once gmsh.initialize() was executed.
    # Use abs-path instead.
    if elmergrid == 'win':
        elmergrid = 'C:/Program Files/Elmer 8.4-Release/bin/ElmerGrid.exe'

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


def run_elmer_solver(sim_dir, elmersolver='win'):
    """Run ElmerSolver with input file case.sif

    Args:
        sim_dir (str): Simulation directory
        elmersolver (str, optional): ElmerSolver executable
    """
    # On Windows ElmerSolver.exe is not found once gmsh.initialize() was executed.
    # Use abs-path instead.
    if elmersolver == 'win':
        elmersolver = 'C:/Program Files/Elmer 8.4-Release/bin/ElmerSolver.exe'

    args = [elmersolver, 'case.sif']
    with open(sim_dir + '/elmersolver.log', 'w') as f:
        subprocess.run(args, cwd=sim_dir, stdout=f, stderr=f)


def scan_logfile(sim_dir):
    """Scan log file for errors and warnings.

    Args:
        sim_dir (str): Simulation directory

    Returns:
        list[str], list[str]: error messages, warnings
    """
    with open(sim_dir + '/elmersolver.log', 'r') as f:
        log = f.readlines()
    print(log)
    err = []
    warn = []
    for line in log:
        if 'ERROR' in line:
            err.append(line[:-1])
        if 'WARNING' in line:
            err.append(line[:-1])
    return err, warn


if __name__ == "__main__":
    err, warn = scan_logfile('./simdata')
    print(err)
    print(warn)
