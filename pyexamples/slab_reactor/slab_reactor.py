# pyexamples/slab_reactor/slab_reactor.py
#
#   Assembly     kinf
#   ========   =========
#      A       1.330098
#      B       1.299289
#      C       0.679513
#      D       0.191268
#
#   Core     keff
#   =====   =========
#     1     1.258874
#     2     1.006969
#     3     0.804291
#

import numpy as np
import time
from detran import *
import slab_reactor_materials

#------------------------------------------------------------------------------#
# Input
#------------------------------------------------------------------------------#

inp = InputDB.Create()
inp.put_int("number_groups",   2)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 1)
inp.put_dbl("inner_tolerance", 1e-4)
inp.put_int("inner_print_out", 0)
inp.put_int("outer_max_iters", 1)
inp.put_dbl("outer_tolerance", 1e-4)
inp.put_int("outer_print_out", 0)
inp.put_int("eigen_max_iters", 1000)
inp.put_dbl("eigen_tolerance", 1e-4)
inp.put_str("bc_left",         "vacuum")
inp.put_str("bc_right",        "vacuum")

#------------------------------------------------------------------------------#
# Material
#------------------------------------------------------------------------------#

mat = slab_reactor_materials.get_materials()

#------------------------------------------------------------------------------#
# Mesh
#------------------------------------------------------------------------------#

# Assembly coarse mesh edges
base = [1.1580, 4.4790, 7.8000, 11.1210, 14.4420, 15.6000]

# and fine mesh counts.
factor = 2
basef = [1, 2, 2, 2, 2, 1]
 
# Several such assemblies to make the total coarse mesh definition
xcm_c = [0.0]
xfm_c = []
for i in range(0, 7) :
  for j in range(0, len(base)) :
    xcm_c.append(base[j] + 15.6*float(i))
    xfm_c.append(basef[j] * factor)

# Make one to do single assemblies, too.
xcm_a  = [0.0000, 1.1580, 4.4790, 7.8000, 11.1210, 14.4420, 15.6000]
xfm_a  = basef;

# Assembly types
assem = [[ 0, 1, 2, 2, 1, 0 ], \
         [ 0, 1, 1, 1, 1, 0 ], \
         [ 0, 1, 3, 3, 1, 0 ], \
         [ 0, 3, 3, 3, 3, 0 ]]

# Cores 1, 2 and 3
core_1 = assem[0]+assem[1]+assem[0]+assem[1]+assem[0]+assem[1]+assem[0]
core_2 = assem[0]+assem[2]+assem[0]+assem[2]+assem[0]+assem[2]+assem[0]
core_3 = assem[0]+assem[3]+assem[0]+assem[3]+assem[0]+assem[3]+assem[0]
core_4 = assem[0]+assem[3]+assem[3]+assem[3]+assem[3]+assem[3]+assem[0]
mesh = Mesh1D.Create(xfm_c, xcm_c, core_3);
mesh.display()
mesh_ref = Mesh1D(xfm_c, xcm_c, core_1);

#------------------------------------------------------------------------------#
# Quadratures
#------------------------------------------------------------------------------#

quad = GaussLegendre.Create(4)
quad.display()

#------------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------------#

# State
state = State.Create(inp, mesh, quad)
# Empty source
q_e = ExternalSourceSP()
# Fission source
q_f = FissionSource.Create(state, mesh, mat)
q_f.initialize()
# boundary
bound = Boundary1D.Create(inp, mesh, quad)

#------------------------------------------------------------------------------#
# Solve
#------------------------------------------------------------------------------#

solver = PowerIteration1D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)
start = time.time()
solver.solve() 
elapsed = (time.time() - start)
print elapsed, " seconds"

#------------------------------------------------------------------------------#
# Post process
#------------------------------------------------------------------------------#

# Plots and things
v1 = np.asarray(state.phi(0))
v2 = np.asarray(state.phi(1))
qf = np.asarray(q_f.density())
#print v
#mesh_ref.plot_mesh_map("MATERIAL")
#mesh_ref.plot_flux(v2)
#mesh_ref.plot_flux(v2)
##print dir(State)


