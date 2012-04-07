# pyexamples/simple_box.py
#
# A simple 2-d square region, 1 group, uniform isotropic source

import detran
print dir(detran)

from detran import *

# Material
mat = Material.Create(1, 1, False)
mat.set_sigma_t(0, 0, 1.0)
mat.set_sigma_s(0, 0, 0, 0.1)
mat.display()

# Mesh
cm = [0.0, 5.0]
fm = [10]
cm_mat = [0]
mesh = Mesh2D.Create(fm, fm, cm, cm, cm_mat)
print mesh.number_cells()
print mesh
# Quadrature
quad = QuadrupleRange.Create(2)

# Input
inp = InputDB.Create()
inp.put_int("number_groups", 1)

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





