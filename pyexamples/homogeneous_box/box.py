# pyexamples/box.py
#
# A simple 2-d square region, 1 or 7 group reactor.  The goal here is
# to provide a very simple benchmark for the Denovo-DGM effort.

import numpy as np
import time
from detran import *

#-----------------------------------------------------------------------------#
# INPUT
#-----------------------------------------------------------------------------#

inp = InputDB.Create()
inp.put_int("number_groups",      7)
inp.put_str("equation",           "sc")
inp.put_int("inner_max_iters",    1)
inp.put_dbl("inner_tolerance",    1e-7)
inp.put_int("inner_print_out",    0)
inp.put_int("outer_max_iters",    1)
inp.put_dbl("outer_tolerance",    1e-7)
inp.put_int("outer_print_out",    0)
inp.put_int("eigen_max_iters",    10000)
inp.put_dbl("eigen_tolerance",    1e-9)
inp.put_str("bc_left",            "reflect")
inp.put_str("bc_right",           "reflect")
inp.put_str("bc_bottom",          "reflect")
inp.put_str("bc_top",             "vacuum")
inp.put_int("store_angular_flux", 1)
inp.display()

#-----------------------------------------------------------------------------#
# MATERIAL
#-----------------------------------------------------------------------------#

if inp.get_int("number_groups") == 1 :

  mat = Material.Create(1, 1, False)
  mat.set_sigma_t(0, 0,    1.0)
  mat.set_sigma_s(0, 0, 0, 0.5)
  mat.set_nu_sigma_f(0, 0, 0.5)
  mat.set_chi(0, 0,        1.0)

else :

  # Material
  mat = Material.Create(7, 1, True)
  m = 0
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.77949E-01)
  mat.set_sigma_t(m, 1, 3.29805E-01)
  mat.set_sigma_t(m, 2, 4.80388E-01)
  mat.set_sigma_t(m, 3, 5.54367E-01)
  mat.set_sigma_t(m, 4, 3.11801E-01)
  mat.set_sigma_t(m, 5, 3.95168E-01)
  mat.set_sigma_t(m, 6, 5.64406E-01)
  # Fission times nu
  mat.set_nu_sigma_f(m, 0, 7.21206E-03*2.78145E+00)
  mat.set_nu_sigma_f(m, 1, 8.19301E-04*2.47443E+00)
  mat.set_nu_sigma_f(m, 2, 6.45320E-03*2.43383E+00)
  mat.set_nu_sigma_f(m, 3, 1.85648E-02*2.43380E+00)
  mat.set_nu_sigma_f(m, 4, 1.78084E-02*2.43380E+00)
  mat.set_nu_sigma_f(m, 5, 8.30348E-02*2.43380E+00)
  mat.set_nu_sigma_f(m, 6, 2.16004E-01*2.43380E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 0 <- g'
  mat.set_sigma_s(m, 0, 0, 1.27537E-01)
  # 1 <- g'
  mat.set_sigma_s(m, 1, 0, 4.23780E-02)
  mat.set_sigma_s(m, 1, 1, 3.24456E-01)
  # 2 <- g'
  mat.set_sigma_s(m, 2, 0, 9.43740E-06)
  mat.set_sigma_s(m, 2, 1, 1.63140E-03)
  mat.set_sigma_s(m, 2, 2, 4.50940E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 3, 0, 5.51630E-09)
  mat.set_sigma_s(m, 3, 1, 3.14270E-09)
  mat.set_sigma_s(m, 3, 2, 2.67920E-03)
  mat.set_sigma_s(m, 3, 3, 4.52565E-01)
  mat.set_sigma_s(m, 3, 4, 1.25250E-04*0.0)
  # 4 <- g'
  mat.set_sigma_s(m, 4, 3, 5.56640E-03)
  mat.set_sigma_s(m, 4, 4, 2.71401E-01)
  mat.set_sigma_s(m, 4, 5, 1.29680E-03*0.0)
  # 5 <- g'
  mat.set_sigma_s(m, 5, 4, 1.02550E-02)
  mat.set_sigma_s(m, 5, 5, 2.65802E-01)
  mat.set_sigma_s(m, 5, 6, 8.54580E-03*0.0)
  # 6 <- g'
  mat.set_sigma_s(m, 6, 4, 1.00210E-08)
  mat.set_sigma_s(m, 6, 5, 1.68090E-02)
  mat.set_sigma_s(m, 6, 6, 2.73080E-01)
  
mat.finalize()
mat.display()

#-----------------------------------------------------------------------------#
# MESH
#-----------------------------------------------------------------------------#

cm        = [0.0, 10.0]
fm        = [5]
cm_mat    = [0]
mesh      = Mesh2D.Create(fm, fm, cm, cm, cm_mat)
mesh_ref  = Mesh2D(fm, fm, cm, cm, cm_mat)

#-----------------------------------------------------------------------------#
# QUADRATURE
#-----------------------------------------------------------------------------#

quad = LevelSymmetric.Create(8, 2)
quad.display()

#-----------------------------------------------------------------------------#
# MISC. SETUP
#-----------------------------------------------------------------------------#

# State
state = State.Create(inp, mesh, quad)
# Empty external source.
q_e = ExternalSourceSP()
# Fission source
q_f = FissionSource.Create(state, mesh, mat)
q_f.initialize()
# Boundary
bound = Boundary2D.Create(inp, mesh, quad)

#-----------------------------------------------------------------------------#
# SOLVE
#-----------------------------------------------------------------------------#

solver = PowerIteration2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)

start = time.time()
solver.solve()
elapsed = (time.time() - start)

print elapsed, " seconds"

phi_0 = np.asarray(state.phi(0))

psi_0= np.asarray(state.psi(0, 0, 0))
psi_1= np.asarray(state.psi(0, 0, 1))
psi_2= np.asarray(state.psi(0, 0, 2))
psi_3= np.asarray(state.psi(0, 0, 3))

for i in range(0, 25):
  print '%10.8f %10.8f %10.8f %10.8f %10.8f' % (phi_0[i]/phi_0[0], psi_0[i]/psi_0[0], psi_1[i]/psi_0[0], psi_2[i]/psi_0[0], psi_3[i]/psi_0[0])

mesh_ref.plot_flux(phi_0)
quad.display()
