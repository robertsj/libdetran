# pyexamples/pin_cell.py
#
# 2-D pin cell
import numpy as np
import time
from detran import *

pin = PinCell.Create(1.27, [0.30,0.40,0.50], [0,0,0,1])
pin.meshify(300)
mesh = pin.mesh()
mesh_ref = pin.mesh_ref()

# Input
inp = InputDB.Create()
inp.put_int("number_groups",   1)
inp.put_str("equation",        "sc")
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

# Quadrature
quad = QuadrupleRange.Create(8)

# State
state = State.Create(inp, pin.mesh(), quad)

# Constant source
q_e = ExternalSourceSP()

# Uninitialized fission source
q_f = FissionSource.Create(state, pin.mesh(), mat)
q_f.initialize()

# boundary
bound = Boundary2D.Create(inp, pin.mesh(), quad)
solver = PowerIteration2D.Create(inp, state, pin.mesh(), mat, quad, bound, q_e, q_f)

# Solve and time.
start = time.time()
solver.solve()
elapsed = (time.time() - start)
print elapsed, " seconds"

v = np.asarray(state.phi(0))
mesh_ref.plot_flux(v)


