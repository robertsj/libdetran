# pyexamples/simple_box.py
#
# A simple slab, uniform isotropic source
import numpy as np
import time
from detran import *
print dir(Mesh1D)
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
inp.put_str("bc_left",         "vacuum")
inp.put_str("bc_right",        "vacuum")

# Material
mat = Material.Create(1, 2, False)
# 0
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.5)
mat.set_nu_sigma_f(0, 0, 0.2)
mat.set_chi(0, 0,        1.0)
# 1
mat.set_sigma_t(1, 0,    1.0)
mat.set_sigma_s(1, 0, 0, 0.5)
mat.set_nu_sigma_f(1, 0, 0.5)
mat.set_chi(1, 0,        1.0)

mat.finalize()
#mat.display()

# Mesh
cm = [0.0, 10.0, 20.0]
fm = [100, 100]
cm_mat = [0, 1]
mesh = Mesh1D.Create(fm, cm, cm_mat)
mesh.display()

mesh_ref = Mesh1D(fm, cm, cm_mat)

# Quadrature
quad = GaussLegendre.Create(32)
quad.display()

# State
state = State.Create(inp, mesh, quad)

# Constant source
q_e = ExternalSourceSP()#ConstantSource.Create(mesh, quad, 1, 1.0)

# Uninitialized fission source
q_f = FissionSource.Create(state, mesh, mat)
q_f.initialize()

# boundary
bound = Boundary1D.Create(inp, mesh, quad)

# solver
#solver = GaussSeidel1D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)
solver = PowerIteration1D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)
start = time.time()
solver.solve() # solve group 0.


elapsed = (time.time() - start)
print elapsed, " seconds"
#v = np.asarray(state.phi(0))
v = np.asarray(q_f.density())
print v

mesh_ref.plot_flux(v)
#print dir(State)


