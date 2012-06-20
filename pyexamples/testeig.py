import numpy as np
import matplotlib.pyplot as plt
from detran import *
import sys
import time
# Input
inp = InputDB.Create()
inp.put_str("problem_type",    "eigenvalue")
inp.put_int("number_groups",   2)
inp.put_str("inner_solver",    "GMRES")
inp.put_int("inner_max_iters", 1000)
inp.put_dbl("inner_tolerance", 1e-9)
inp.put_int("inner_print_out", 0)
inp.put_str("outer_solver",    "GS")
inp.put_int("outer_max_iters", 100)
inp.put_dbl("outer_tolerance", 1e-5)
inp.put_int("outer_print_out", 0)
inp.put_str("eigen_solver",    "SLEPc")
inp.put_dbl("eigen_tolerance", 1e-5)
inp.put_int("eigen_max_iters", 100)
inp.put_str("bc_left",         "vacuum")
inp.put_str("bc_right",        "reflect")
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
fm = [1000]
cm = [0, 50]
mesh = Mesh1D.Create(fm, cm, [0])
# Execute
execute = Execute1D(sys.argv)
execute.initialize(inp, mat, mesh)
execute.solve()
# Postprocess
state = execute.get_state()
#print np.asarray(state.phi(0)) 
#print np.asarray(state.phi(1)) 
p0 = np.asarray(state.phi(0))
p0 = p0 / sum(p0)
p1 = np.asarray(state.phi(1))
p1 = p1 / sum(p1)
#plt.plot(range(0,1000), p0, range(0,1000), p1, linewidth=3)
#plt.grid(True)
#plt.show()
#[ 1.41157171  1.77647353  2.09314514  2.34491177  2.51949292  2.60880394]
#[ 0.1405278   0.24012095  0.31761254  0.37495543  0.41297769  0.43195761]
execute.finalize()
