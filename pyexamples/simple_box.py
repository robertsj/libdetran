# pyexamples/simple_box.py
#
# A simple 2-d square region, 1 group, uniform isotropic source
import numpy as np
import detran
print dir(detran.detran_transport)

from detran import *

# Input
inp = InputDB.Create()
inp.put_int("number_groups",   1)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 20)
inp.put_dbl("inner_tolerance", 1e-8)

# Material
mat = Material.Create(1, 2, False)
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.2)
mat.set_sigma_t(1, 0,    1.0)
mat.set_sigma_s(1, 0, 0, 0.9)
mat.finalize()
mat.display()

# Mesh
#cm = [0.0, 1.0]
#fm = [3]
#cm_mat = [0]
cm = [0.0, 5.0, 10.0]
fm = [  100, 100]
cm_mat = [1, 0, 0, 0]

mesh = Mesh2D.Create(fm, fm, cm, cm, cm_mat)
mesh.display()

mesh2 = Mesh2D(fm, fm, cm, cm, cm_mat)

# Quadrature
quad = QuadrupleRange.Create(50)
quad.display()

# State
state = State.Create(inp, mesh, quad)

# Constant source
q_e = ConstantSource.Create(mesh, quad, 1, 1.0)

# Uninitialized fission source
q_f = FissionSourceSP()

# boundary
bound = Boundary2D.Create(inp, mesh, quad)
print q_e
solver = SourceIteration2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)
solver.solve(0) # solve group 0.

v = np.asarray(state.phi(0))
#print v.reshape(5, 5)
mesh2.plot_flux(v)
print dir(State)


