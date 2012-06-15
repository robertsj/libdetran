import numpy as np
import matplotlib.pyplot as plt
from detran import *
import sys
import time
# Input
inp = InputDB.Create()
inp.put_str("problem_type",    "eigenvalue")
inp.put_str("inner_solver",    "GMRES")
inp.put_int("inner_max_iters", 1000)
inp.put_dbl("inner_tolerance", 1e-9)
inp.put_int("inner_print_out", 1)
inp.put_str("outer_solver",    "GS")
inp.put_int("outer_max_iters", 100)
inp.put_dbl("outer_tolerance", 1e-8)
inp.put_int("outer_print_out", 0)
# Material
mat = Material.Create(2, 1, False)
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.5)
mat.set_sigma_f(0, 0,    0.5)
mat.set_chi(0, 0,        1.0)
mat.finalize()
# Mesh
fm = [6]
cm = [0, 5]
mesh = Mesh1D.Create(fm, cm, [0])
# Execute
execute = Execute1D(sys.argv)
print sys.argv
execute.initialize(inp, mat, mesh)
execute.solve()
# Postprocess
state = execute.get_state()
print np.asarray(state.phi(0)) 
plt.plot(np.asarray(state.phi(0)), linewidth=3)
plt.grid(True)
#plt.show()
#4.75311163  7.95478996  9.00361562  9.00361562  7.95478996  4.75311163
execute.finalize()
