import os
import gmsh
from pyelmer import elmer, post
from pyelmer import execute
from objectgmsh import add_physical_group
from math import floor, log10


############################################################################
### Settings
############################################################################

# The following switches defined are for debugging purposes
# run Elmer solver. Default: True
run_solver = True
# run Gmsh and create mesh. Default: True
run_mesher = True
# run Gmsh GUI in order to view and verify the mesh manually before Elmer simulation. Default: False
run_gmsh_gui = False

# Parameter Definition
w = 6e-3  # pad width
l = 6e-3  # pad length
r = 2.95e-3  # pad fillet radius
t = 0.035e-3  # metal thickness
tol = 1e-6  # search tol
hsub = 1.524e-3  # RO4003C substrate height
hmb = 0.1e-3  # bottom metal height
wsub = 4 * max([w, l])  # substrate width
h_ab = 5 * (hsub + hmb + t)  # air box height
w_ab = 1.2 * wsub  # air box width

# set up working subdirectory
sim_dir = "./simdata/"

if not os.path.exists(sim_dir):
    os.mkdir(sim_dir)

############################################################################
### Some useful functions
############################################################################


def getMyEntitiesInBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax, dim):
    entities = gmsh.model.getEntitiesInBoundingBox(
        xmin, ymin, zmin, xmax, ymax, zmax, dim
    )
    all_entities = []
    for item in entities:
        all_entities.append(item[1])
    return all_entities


# scan err, war, stats and results
def extract_results_logfile(sim_dir):
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
    for line in log:
        if "StatElecSolve:  Capacitance" in line:  # extract capacitance from log file
            s = " ".join(line.split()).split(" ")
            capacitance = float(s[3])
        if (
            "StatElecSolve:  Relative Change" in line
        ):  # extract relative change from log file
            s = " ".join(line.split()).split(" ")
            rel_change = float(s[4])
    return capacitance, rel_change


# engineering format
def powerise10(x):
    """Returns x as a*10**b with 0 <= a < 10"""
    if x == 0:
        return 0, 0
    Neg = x < 0
    if Neg:
        x = -x
    a = 1.0 * x / 10 ** (floor(log10(x)))
    b = int(floor(log10(x)))
    if Neg:
        a = -a
    return a, b


def eng(x):
    """Return a string representing x in an engineer friendly notation"""
    a, b = powerise10(x)
    if -3 < b < 3:
        return "%.4g" % x
    a = a * 10 ** (b % 3)
    b = b - b % 3
    return "%.4gE%s" % (a, b)


############################################################################
### Geometry modeling using gmsh
############################################################################

gmsh.initialize()

gmsh.model.add("3d_electrostatic_capacitance")
geom = gmsh.model.occ

# top metal
m1 = geom.addBox(-w / 2, -l / 2, 0, w, l, t)
geom.fillet([m1], [1, 3, 5, 7], [r], True)

# substrate
sub = geom.addBox(-wsub / 2, -wsub / 2, -hsub, wsub, wsub, hsub)

# bottom metal
m2 = geom.addBox(-wsub / 2, -wsub / 2, -(hsub + hmb), wsub, wsub, hmb)

# airbox: we mesh only air and substrate, but not metal
# design proper height of air box to cover fringing fields
ab = geom.addBox(-w_ab / 2, -w_ab / 2, -(hsub + 2 * hmb), w_ab, w_ab, h_ab)

# remove metal volumes
geom.cut([(3, ab)], [(3, m1), (3, m2)], removeObject=True, removeTool=True)

geom.synchronize()
geom.fragment([(3, sub)], [(3, ab)])

##########################################################
# structured meshing would be advantageous but needs further optimization

# ## M1
# NN = 4
# tf_lines  = [2,13,21,23,24,22,15,3]
# for k in tf_lines:
#     gmsh.model.mesh.setTransfiniteCurve(k, NN)

# ## Sub
# NN = int((w-2*r)*10*1e3)
# tf_lines = [7,16,8,17]
# for k in tf_lines:
#     gmsh.model.mesh.setTransfiniteCurve(k, NN)

# ## MB
# NN = int((l-2*r)*10*1e3)
# tf_lines = [1,4,11,20]
# for k in tf_lines:
#     gmsh.model.mesh.setTransfiniteCurve(k, NN)

# NN = 5
# tf_lines = [42,45,47,41,50,53,55,49]
# for k in tf_lines:
#     gmsh.model.mesh.setTransfiniteCurve(k, NN)

# NN = int(wsub*50*1e3)
# tf_lines = [38,44,52,37,43,51,39,46,54,40,48,56]
# for k in tf_lines:
#     gmsh.model.mesh.setTransfiniteCurve(k, NN)
##########################################################

############################################################################
### Physical Groups and Boundary Conditions
############################################################################
geom.synchronize()

# these two were identified manually in Gmsh GUI
# volumes = gmsh.model.getEntities(dim=3) # check volume numbers after fragment
ph_sub = add_physical_group(3, [2], "substrate")
ph_ab = add_physical_group(3, [3], "airbox")

# the others are identified using own function 'getMyEntitiesInBoundingBox
m1_sfs = getMyEntitiesInBoundingBox(
    -w / 2 - tol, -l / 2 - tol, 0 - tol, w / 2 + tol, l / 2 + tol, t + tol, 2
)
ph_m1_sfs = add_physical_group(2, m1_sfs, "metal1")

m2_sfs = getMyEntitiesInBoundingBox(
    -wsub / 2 - tol,
    -wsub / 2 - tol,
    -(hsub + hmb) - tol,
    wsub / 2 + tol,
    wsub / 2 + tol,
    -hsub + tol,
    2,
)
ph_m2_sfs = add_physical_group(2, m2_sfs, "metal2")

sub_sfs = getMyEntitiesInBoundingBox(
    -wsub / 2 - tol,
    -wsub / 2 - tol,
    -hsub - tol,
    wsub / 2 + tol,
    wsub / 2 + tol,
    tol,
    2,
)

temp_sfs = getMyEntitiesInBoundingBox(
    -w_ab / 2 - tol,
    -w_ab / 2 - tol,
    -h_ab / 2 - tol,
    w_ab / 2 + tol,
    w_ab / 2 + tol,
    h_ab / 2 + tol,
    2,
)
ab_sfs = [x for x in temp_sfs if x not in (m1_sfs + sub_sfs + m2_sfs)]
ph_ab_sfs = add_physical_group(2, ab_sfs, "ab_sfs")

############################################################################
### Meshing
############################################################################

# We can activate the calculation of mesh element sizes based on curvature
# (here with a target of 90 elements per 2*Pi radians):
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 90)

# Finally we apply an elliptic smoother to the grid to have a more regular
# mesh:
gmsh.option.setNumber("Mesh.Smoothing", 10)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # faster
# gmsh.option.setNumber('General.NumThreads', 8)
# gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2e-3)

if run_mesher:
    geom.synchronize()
    gmsh.model.mesh.generate(dim=3)
    gmsh.write(sim_dir + "/3d_electrostatic_capacitance.msh")
    # Preview mesh.
if run_gmsh_gui:
    gmsh.fltk.run()

# Clear mesh and close gmsh API.
gmsh.clear()
gmsh.finalize()

############################################################################
### Elmer Setup
############################################################################

sim = elmer.load_simulation("3D_steady", "my_simulations.yml")
# adding constants is very important, otherwise the solver calculates wrong results!
sim.constants.update({"Permittivity of Vacuum": "8.8542e-12"})
sim.constants.update({"Gravity(4)": "0 -1 0 9.82"})
sim.constants.update({"Boltzmann Constant": "1.3807e-23"})
sim.constants.update({"Unit Charge": "1.602e-19"})

# materials
air = elmer.load_material("air", sim, "my_materials.yml")
ro4003c = elmer.load_material("ro4003c", sim, "my_materials.yml")

# solver
solver_electrostatic = elmer.load_solver("Electrostatics", sim, r"my_solvers.yml")
# very important, the value must match the boundary condition abs(potential difference) !!!
# otherwise the capacitance will be calculated wrong !
solver_electrostatic.data.update({"Potential Difference": "1.0"})

# equation
eqn = elmer.Equation(sim, "main", [solver_electrostatic])

# bodies
bdy_sub = elmer.Body(sim, "substrate", [ph_sub])
bdy_sub.material = ro4003c
bdy_sub.equation = eqn

bdy_ab = elmer.Body(sim, "airbox", [ph_ab])
bdy_ab.material = air
bdy_ab.equation = eqn

# boundaries
bndry_m1 = elmer.Boundary(sim, "top metal", [ph_m1_sfs])
bndry_m1.data.update({"Potential": "1.0"})

bndry_m2 = elmer.Boundary(sim, "bottom metal", [ph_m2_sfs])
bndry_m2.data.update({"Potential": "0.0"})

bndry_airbox = elmer.Boundary(sim, "FarField", [ph_ab_sfs])
bndry_airbox.data.update({"Electric Infinity BC": "True"})

# export
sim.write_startinfo(sim_dir)
sim.write_sif(sim_dir)

if run_mesher:
    execute.run_elmer_grid(sim_dir, "3d_electrostatic_capacitance.msh")

##############
# execute ElmerGrid & ElmerSolver
if run_solver:
    execute.run_elmer_solver(sim_dir)
    ###############
    # scan log for errors and warnings
    err, warn, stats = post.scan_logfile(sim_dir)
    capacitance, rel_change = extract_results_logfile(sim_dir)
    print("## RESULTS BEGIN #############################")
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)
    print("Relative Change:", f"{rel_change:.2E}")
    print("##############################################")
    print("w:", eng(w))
    print("l:", eng(l))
    print("r:", eng(r))
    print("Capacitance:", eng(capacitance))
    print("##############################################")
