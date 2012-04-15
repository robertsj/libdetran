# pyexamples/simple_box.py
#
# A simple 2-d square region, 1 group, uniform isotropic source
import numpy as np
import time


import detran
print dir(detran.detran_solvers)

from detran import *

# Input
inp = InputDB.Create()
inp.put_int("number_groups",   1)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 1)
inp.put_dbl("inner_tolerance", 1e-4)
inp.put_int("inner_print_out", 0)
inp.put_int("outer_max_iters", 1)
inp.put_dbl("outer_tolerance", 1e-4)
inp.put_int("outer_print_out", 0)
inp.put_int("eigen_max_iters", 1000)
inp.put_dbl("eigen_tolerance", 1e-4)
inp.put_str("bc_left",         "reflect")
inp.put_str("bc_right",        "reflect")
inp.put_str("bc_bottom",       "reflect")
inp.put_str("bc_top",          "reflect")


# Material
mat = Material.Create(1, 2, False)
#
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.5)
mat.set_nu_sigma_f(0, 0, 0.5)
mat.set_chi(0, 0,        1.0)
#
mat.set_sigma_t(1, 0,    1.0)
mat.set_sigma_s(1, 0, 0, 0.9)
mat.finalize()
#mat.display()

# Mesh
#cm = [0.0, 10.0]
#fm = [100]
#cm_mat = [0]
cm = [0.0, 5.0, 10.0, 15.0, 20.0]
fm = [ 20, 20, 20, 20]
cm_mat = [1, 1, 1, 0,
          0, 0, 0, 0,
          1, 0, 0, 0,
          0, 0, 0, 0]

cm_mat2 = [[1, 1, 1, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]]
print cm_mat2
print np.flipud(cm_mat2)
print np.transpose(np.flipud(cm_mat2))
print np.squeeze(np.asarray(np.reshape(np.transpose(np.flipud(cm_mat2)), (16,-1))))

mesh = Mesh2D.Create(fm, fm, cm, cm, cm_mat)
mesh.display()

mesh2 = Mesh2D(fm, fm, cm, cm, cm_mat)

# Quadrature
quad = QuadrupleRange.Create(8)
#quad.display()

# State
state = State.Create(inp, mesh, quad)

# Constant source
q_e = ExternalSourceSP()#ConstantSource.Create(mesh, quad, 1, 1.0)

# Uninitialized fission source
q_f = FissionSource.Create(state, mesh, mat)
q_f.initialize()

# boundary
bound = Boundary2D.Create(inp, mesh, quad)
print(dir(bound))
bv = bound(0, 0, 0, 0)
print bv
solver = PowerIteration2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)

start = time.time()
solver.solve() # solve group 0.
elapsed = (time.time() - start)
print elapsed, " seconds"

v = np.asarray(state.phi(0))
#print v[0:3]
pin2.plot_flux(v)


#print dir(State)


