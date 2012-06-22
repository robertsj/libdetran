import numpy as np
import time
import sys
from detran import *
#
from assemblies_c5g7 import get_assemblies
from pins_c5g7 import get_pins
from material_c5g7 import get_materials
from plot_utils import *
#-----------------------------------------------------------------------------#
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_str("equation",                 "dd")
inp.put_str("problem_type",             "eigenvalue")
#
inp.put_str("inner_solver",             "GMRES")
inp.put_int("inner_max_iters",          1)
inp.put_dbl("inner_tolerance",          1e-8)
inp.put_int("inner_print_out",          1)
inp.put_int("inner_print_interval",     10)
#
inp.put_str("outer_solver",             "KrylovMG")
inp.put_int("outer_max_iters",          0)
inp.put_dbl("outer_tolerance",          1e-8)
inp.put_int("outer_print_out",          1)
inp.put_int("outer_print_interval",     2)
#inp.put_int("outer_upscatter_cutoff",   0) 
# 0.547002225
inp.put_str("eigen_solver",             "SLEPc")
inp.put_int("eigen_max_iters",          1000)
inp.put_dbl("eigen_tolerance",          1e-14)
inp.put_int("eigen_print_out",          2)
inp.put_int("eigen_print_interval",     1)
#
inp.put_str("bc_left",                  "vacuum")
inp.put_str("bc_right",                 "vacuum")
inp.put_str("bc_bottom",                "vacuum")
inp.put_str("bc_top",                   "vacuum")
#
inp.put_str("quad_type",                "quadruplerange")
inp.put_int("quad_order",               2)
#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
mat = get_materials()

#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
assemblies = get_assemblies(7, True)
mesh = assemblies[0].mesh()
#-----------------------------------------------------------------------------#
# Execute
#-----------------------------------------------------------------------------#
execute = Execute2D(sys.argv)
execute.initialize(inp, mat, mesh)
t = time.time()
execute.solve()
print "elapsed = ", time.time()-t
#-----------------------------------------------------------------------------#
# Post-Process
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
execute.finalize()

