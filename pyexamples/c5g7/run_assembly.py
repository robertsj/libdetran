import numpy as np
import time
from detran import *
import assemblies_c5g7
import material_c5g7

assemblies= assemblies_c5g7.get_assemblies()
assem1 = assemblies[0]
mesh = assem1.mesh()
meshref = assem1.mesh_ref()
matmap = mesh.mesh_map("MATERIAL")
#meshref.plot_mesh_map("MATERIAL")

#pins = assemblies_c5g7.get_pins()
#pin  = pins[0]
#mesh = pin.mesh()

#cm = [0.0, 1.0]
#fm = [  10]
#cm_mat = [6]
#mesh = Mesh2D.Create(fm, fm, cm, cm, cm_mat)
#meshref = Mesh2D(fm, fm, cm, cm, cm_mat)
# Input
inp = InputDB.Create()
inp.put_int("number_groups",   7)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 1)
inp.put_dbl("inner_tolerance", 1e-9)
inp.put_int("inner_print_out", 0)
inp.put_int("outer_max_iters", 1)
inp.put_dbl("outer_tolerance", 1e-9)
inp.put_int("outer_print_out", 0)
inp.put_int("eigen_max_iters", 10)
inp.put_dbl("eigen_tolerance", 1e-6)
inp.put_str("bc_left",         "reflect")
inp.put_str("bc_right",        "reflect")
inp.put_str("bc_bottom",       "reflect")
inp.put_str("bc_top",          "reflect")

# Material
mat = material_c5g7.get_materials()

# Quadrature
quad = QuadrupleRange.Create(8)

# State
state = State.Create(inp, mesh, quad)

# Constant source
q_e = ExternalSourceSP()
#ConstantSource.Create(mesh, quad, 7, 1.0)

# Uninitialized fission source
q_f = FissionSource.Create(state, mesh, mat)
q_f.initialize()

# boundary
bound = Boundary2D.Create(inp, mesh, quad)
solver = PowerIteration2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)
#solver = GaussSeidel2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)

# Solve and time.
start = time.time()
solver.solve()
elapsed = (time.time() - start)
print elapsed, " seconds"

#vv = np.zeros(7)
#for i in range(0, 7):
#  v = np.asarray(state.phi(i))
#  vv[i] = v[0]

#print vv
v = np.asarray(q_f.density())
#meshref = pin.mesh_ref()
print v
meshref.plot_flux(v)
