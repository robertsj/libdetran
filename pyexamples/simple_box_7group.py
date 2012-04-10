# pyexamples/simple_box.py
#
# A simple 2-d square region, 1 group, uniform isotropic source
import numpy as np
import time


import detran
print dir(detran.detran_solvers)

from detran import *

# Input
inp = InputDB.Create()
inp.put_int("number_groups",   7)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 10)
inp.put_dbl("inner_tolerance", 1e-4)
inp.put_int("outer_max_iters", 3)
inp.put_dbl("outer_tolerance", 1e-4)

inp.put_str("bc_left",         "reflect")
inp.put_str("bc_right",        "vacuum")
inp.put_str("bc_bottom",       "vacuum")
inp.put_str("bc_top",          "vacuum")


# Material
mat = Material.Create(7, 2, False)
# --------------------------------------------
# Material 0: UO2 fuel-clad
# --------------------------------------------
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
# 1 <- g'
mat.set_sigma_s(m, 0, 0, 1.27537E-01)
# 2 <- g'
mat.set_sigma_s(m, 1, 0, 4.23780E-02)
mat.set_sigma_s(m, 1, 1, 3.24456E-01)
# 3 <- g'
mat.set_sigma_s(m, 2, 0, 9.43740E-06)
mat.set_sigma_s(m, 2, 1, 1.63140E-03)
mat.set_sigma_s(m, 2, 2, 4.50940E-01)
# 4 <- g'
mat.set_sigma_s(m, 3, 0, 5.51630E-09)
mat.set_sigma_s(m, 3, 1, 3.14270E-09)
mat.set_sigma_s(m, 3, 2, 2.67920E-03)
mat.set_sigma_s(m, 3, 3, 4.52565E-01)
mat.set_sigma_s(m, 3, 4, 1.25250E-04)
# 5 <- g'
mat.set_sigma_s(m, 4, 3, 5.56640E-03)
mat.set_sigma_s(m, 4, 4, 2.71401E-01)
mat.set_sigma_s(m, 4, 5, 1.29680E-03)
# 6 <- g'
mat.set_sigma_s(m, 5, 4, 1.02550E-02)
mat.set_sigma_s(m, 5, 5, 2.65802E-01)
mat.set_sigma_s(m, 5, 6, 8.54580E-03)
# 7 <- g'
mat.set_sigma_s(m, 6, 4, 1.00210E-08)
mat.set_sigma_s(m, 6, 5, 1.68090E-02)
mat.set_sigma_s(m, 6, 6, 2.73080E-01)
# --------------------------------------------
# Material 1: Moderator
# --------------------------------------------
m = 1
# Transport cross section
mat.set_sigma_t(m, 0, 1.59206E-01)
mat.set_sigma_t(m, 1, 4.12970E-01)
mat.set_sigma_t(m, 2, 5.90310E-01)
mat.set_sigma_t(m, 3, 5.84350E-01)
mat.set_sigma_t(m, 4, 7.18000E-01)
mat.set_sigma_t(m, 5, 1.25445E+00)
mat.set_sigma_t(m, 6, 2.65038E+00)
# Scattering
# 1 <- g'
mat.set_sigma_s(m, 0, 0, 4.44777E-02)
# 2 <- g
mat.set_sigma_s(m, 1, 0, 1.13400E-01)
mat.set_sigma_s(m, 1, 1, 2.82334E-01)
# 3 <- g'
mat.set_sigma_s(m, 2, 0, 7.23470E-04)
mat.set_sigma_s(m, 2, 1, 1.29940E-01)
mat.set_sigma_s(m, 2, 2, 3.45256E-01)
# 4 <- g'
mat.set_sigma_s(m, 3, 0, 3.74990E-06)
mat.set_sigma_s(m, 3, 1, 6.23400E-04)
mat.set_sigma_s(m, 3, 2, 2.24570E-01)
mat.set_sigma_s(m, 3, 3, 9.10284E-02)
mat.set_sigma_s(m, 3, 4, 7.14370E-05)
# 5 <- g'
mat.set_sigma_s(m, 4, 0, 5.31840E-08)
mat.set_sigma_s(m, 4, 1, 4.80020E-05)
mat.set_sigma_s(m, 4, 2, 1.69990E-02)
mat.set_sigma_s(m, 4, 3, 4.15510E-01)
mat.set_sigma_s(m, 4, 4, 1.39138E-01)
mat.set_sigma_s(m, 4, 5, 2.21570E-03)
# 6 <- g'
mat.set_sigma_s(m, 5, 1, 7.44860E-06)
mat.set_sigma_s(m, 5, 2, 2.64430E-03)
mat.set_sigma_s(m, 5, 3, 6.37320E-02)
mat.set_sigma_s(m, 5, 4, 5.11820E-01)
mat.set_sigma_s(m, 5, 5, 6.99913E-01)
mat.set_sigma_s(m, 5, 6, 1.32440E-01)
# 7 <- g'
mat.set_sigma_s(m, 6, 1, 1.04550E-06)
mat.set_sigma_s(m, 6, 2, 5.03440E-04)
mat.set_sigma_s(m, 6, 3, 1.21390E-02)
mat.set_sigma_s(m, 6, 4, 6.12290E-02)
mat.set_sigma_s(m, 6, 5, 5.37320E-01)
mat.set_sigma_s(m, 6, 6, 2.48070E+00)
mat.finalize()

# Mesh
#cm = [0.0, 1.0]
#fm = [3]
#cm_mat = [0]
cm = [0.0, 5.0, 10.0, 15.0, 20.0]
fm = [ 10, 10, 10, 10]
cm_mat = [1, 1, 1, 1,
          1, 1, 1, 1,
          1, 1, 1, 1,
          1, 1, 1, 1]

#cm_mat2 = [[1, 1, 1, 0],
#           [0, 0, 0, 0],
#           [0, 0, 0, 0],
#           [0, 0, 0, 0]]
#print cm_mat2
#print np.flipud(cm_mat2)
#print np.transpose(np.flipud(cm_mat2))
#print np.squeeze(np.asarray(np.reshape(np.transpose(np.flipud(cm_mat2)), (16,-1))))

mesh = Mesh2D.Create(fm, fm, cm, cm, cm_mat)
mesh.display()

mesh2 = Mesh2D(fm, fm, cm, cm, cm_mat)

# Quadrature
quad = QuadrupleRange.Create(50)
#quad.display()

# State
state = State.Create(inp, mesh, quad)

# Constant source
q_e = ConstantSource.Create(mesh, quad, 7, 1.0)

# Uninitialized fission source
q_f = FissionSourceSP()

# boundary
bound = Boundary2D.Create(inp, mesh, quad)
print(dir(bound))
bv = bound(0, 0, 0, 0)
print bv
solver = GaussSeidel2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)

start = time.time()
solver.solve() # solve group 0.
elapsed = (time.time() - start)
print elapsed, " seconds"

for i in range(0, 7):
  v = np.asarray(state.phi(i))
  print v[0:3]
  mesh2.plot_flux(v)
#print dir(State)


