import numpy as np
import matplotlib.pyplot as plt
from detran import *
import sys
import time
# Input
inp = InputDB.Create()
inp.put_str("problem_type",             "fixed")
inp.put_str("equation",                 "sd")
#
inp.put_str("inner_solver",             "GMRES")
inp.put_int("inner_max_iters",          1)
inp.put_dbl("inner_tolerance",          1e-6)
inp.put_int("inner_print_out",          1)
inp.put_int("inner_print_interval",     10)
#
inp.put_str("outer_solver",             "KrylovMG")
inp.put_int("outer_max_iters",          2000)
inp.put_dbl("outer_tolerance",          1e-8)
inp.put_int("outer_print_out",          2)
inp.put_int("outer_print_interval",     2)
inp.put_int("outer_upscatter_cutoff",   0)
#
inp.put_str("bc_left",                  "reflect")
inp.put_str("bc_right",                 "reflect")
#
inp.put_int("quad_order",               2)
# Material
mat = Material.Create(2, 1, False)
# Total
mat.set_sigma_t(0, 0, 0.1890);       # (obj, matid, g, value);
mat.set_sigma_t(0, 1, 1.4633);       
# Scattering
f=1
mat.set_sigma_s(0, 0, 0, f*0.1507);    # 1 <- 1
mat.set_sigma_s(0, 0, 1, f*0.0096);    # 1 <- 2
mat.set_sigma_s(0, 1, 0, f*0.0380);    # 2 <- 1
mat.set_sigma_s(0, 1, 1, f*1.4536);    # 2 <- 2
mat.finalize()
#[ 3.47898774  3.97194684  4.12320542  3.97194684  3.47898774]
#[ 0.71103478  0.77493176  0.78600511  0.77493176  0.71103478]
# Mesh
fm = [5]
cm = [0, 10]
mesh = Mesh1D.Create(fm, cm, [0])
execute = Execute1D(sys.argv)
execute.initialize(inp, mat, mesh)
execute.set_external_source(ConstantSource.Create(mesh, execute.get_quadrature(), 2, 1.0))
execute.solve()
# Postprocess
state = execute.get_state()
p0 = np.asarray(state.phi(0))
p1 = np.asarray(state.phi(1))
print p0[0:6]
print p1[0:6]
x  = np.array(range(0, len(p0))) + mesh.dx(0)/2
plt.plot(x, p0, 'k', x, p1, 'b', linewidth=3)
plt.grid(True)
#plt.show()
execute.finalize()
