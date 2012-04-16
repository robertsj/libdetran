# pyexamples/assembly.py
#
# 2-D pin cell
import numpy as np
import time
from detran import *

# Pins
pin1 = PinCell.Create(1.27, [0.45,0.50], [1,1,1])
pin2 = PinCell.Create(1.27, [0.45,0.50], [1,1,1])
pin1.meshify(13)
pin2.meshify(13)

# Assembly
assem = Assembly.Create(3)
assem.add_pincell(pin1)
assem.add_pincell(pin2)
assem.finalize([1,1,1,\
                1,1,1,\
                1,1,1])

#assmesh = assem.mesh_ref()
#matmap = assmesh.mesh_map("MATERIAL")
#assmesh.plot_mesh_map("MATERIAL")

# Input
inp = InputDB.Create()
inp.put_int("number_groups",   2)
inp.put_str("equation",        "sc")
inp.put_int("inner_max_iters", 1)
inp.put_dbl("inner_tolerance", 1e-4)
inp.put_int("inner_print_out", 0)
inp.put_int("outer_max_iters", 1)
inp.put_dbl("outer_tolerance", 1e-4)
inp.put_int("outer_print_out", 0)
inp.put_int("eigen_max_iters", 1000)
inp.put_dbl("eigen_tolerance", 1e-8)
inp.put_str("bc_left",         "reflect")
inp.put_str("bc_right",        "reflect")
inp.put_str("bc_bottom",       "reflect")
inp.put_str("bc_top",          "reflect")

# Material
mat = Material.Create(2, 2, False)

# ---------------------------
# Material 0: Water
# ---------------------------
# Total
mat.set_sigma_t(0, 0, 0.1890) # (obj, matid, g, value)
mat.set_sigma_t(0, 1, 1.4633)
# Fission
mat.set_nu_sigma_f(0, 0, 0.0) # Note, default is zero
mat.set_nu_sigma_f(0, 1, 0.0)
mat.set_chi(0, 0, 0.0)
mat.set_chi(0, 1, 0.0)
# Scattering
mat.set_sigma_s(0, 0, 0, 0.1507) # 1 <- 1
mat.set_sigma_s(0, 0, 1, 0.0000) # 1 <- 2
mat.set_sigma_s(0, 1, 0, 0.0380) # 2 <- 1
mat.set_sigma_s(0, 1, 1, 1.4536) # 2 <- 2
# ---------------------------
# Material 1: Fuel I
# ---------------------------
# Total
mat.set_sigma_t(1, 0, 0.2263) # kinf = 
mat.set_sigma_t(1, 1, 1.0119)
# Fission
mat.set_nu_sigma_f(1, 0, 0.0067)
mat.set_nu_sigma_f(1, 1, 0.1241)
mat.set_chi(1, 0, 1.0)
mat.set_chi(1, 1, 0.0)
# Scattering
mat.set_sigma_s(1, 0, 0, 0.2006) # 1 <- 1
mat.set_sigma_s(1, 0, 1, 0.0000) # 1 <- 2
mat.set_sigma_s(1, 1, 0, 0.0161) # 2 <- 1
mat.set_sigma_s(1, 1, 1, 0.9355) # 2 <- 2

#
mat.finalize()

# Quadrature
quad = QuadrupleRange.Create(50)

# State
state = State.Create(inp, assem.mesh(), quad)

# Constant source
q_e = ExternalSourceSP()

# Uninitialized fission source
q_f = FissionSource.Create(state, assem.mesh(), mat)
q_f.initialize()

# boundary
bound = Boundary2D.Create(inp, assem.mesh(), quad)
solver = PowerIteration2D.Create(inp, state, assem.mesh(), mat, quad, bound, q_e, q_f)

# Solve and time.
start = time.time()
solver.solve()
elapsed = (time.time() - start)
print elapsed, " seconds"

v = np.asarray(state.phi(0))
assem_mesh = assem.mesh_ref()
assem_mesh.plot_flux(v)


