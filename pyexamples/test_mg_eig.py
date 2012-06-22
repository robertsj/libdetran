import numpy as np
import matplotlib.pyplot as plt
from detran import *
import sys
import time
# Input
inp = InputDB.Create()
inp.put_str("problem_type",             "eigenvalue")
inp.put_str("equation",                 "sd")
#
inp.put_str("inner_solver",             "SI")
inp.put_int("inner_max_iters",          1)
inp.put_dbl("inner_tolerance",          1e-6)
inp.put_int("inner_print_out",          1)
inp.put_int("inner_print_interval",     10)
#
inp.put_str("outer_solver",             "GS")
inp.put_int("outer_max_iters",          1000)
inp.put_dbl("outer_tolerance",          1e-8)
inp.put_int("outer_print_out",          2)
inp.put_int("outer_print_interval",     2)
#inp.put_int("outer_upscatter_cutoff",   0)
#
inp.put_str("eigen_solver",             "PI")
inp.put_dbl("eigen_tolerance",          1e-10)
inp.put_int("eigen_max_iters",          1000)
#
inp.put_str("bc_left",                  "vacuum")
inp.put_str("bc_right",                 "reflect")
#
inp.put_int("quad_order",               2)
# Material
mat = Material.Create(2, 1, False)
# Total
mat.set_sigma_t(0, 0, 0.2263);       # (obj, matid, g, value);
mat.set_sigma_t(0, 1, 1.0119);       
# Fission 
mat.set_sigma_f(0, 0, 0.0067);
mat.set_sigma_f(0, 1, 0.1241);   
mat.set_chi(0, 0, 1.0); 
mat.set_chi(0, 1, 0.0);        
# Scattering
mat.set_sigma_s(0, 0, 0, 0.2006);    # 1 <- 1
mat.set_sigma_s(0, 0, 1, 0.0000);    # 1 <- 2
mat.set_sigma_s(0, 1, 0, 0.0161);    # 2 <- 1
mat.set_sigma_s(0, 1, 1, 0.9355);    # 2 <- 2
mat.finalize()
# Mesh
fm = [100]
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
print p0[0:4]
print p1[0:4]
x  = np.array(range(0, len(p0))) + mesh.dx(0)/2
plt.plot(x, p0, 'k', x, p1, 'b', linewidth=3)
plt.grid(True)
#plt.show()
execute.finalize()
