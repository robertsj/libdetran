# pyexamples/performance_tests/simple_box.py
#
# A simple 2-d square region, 1 group, uniform isotropic source.
# Testing OpenMP threading.
# Testing over a range of the number of angles and thread counts.

#-----------------------------------------------------------------------------#
# RUN TEST
#-----------------------------------------------------------------------------#

# Max quad orders
nquad = 5
# How many threads
nproc = 4
# Number of trials
ntime = 5

import numpy as np
import time
import os
import matplotlib.pyplot as plt
import detran

# average times
times = np.zeros(nquad)
# std times
stds  = np.zeros(nquad)
# num angles
nang  = np.zeros(nquad)

for nc in range(0, 1) :

  os.putenv("OMP_NUM_THREADS", str(1))
  os.system("echo $OMP_NUM_THREADS")
  reload(detran)

  #-----------------------------------------------------------------------------#
  # INPUT
  #-----------------------------------------------------------------------------#
  inp = detran.InputDB.Create()
  inp.put_int("number_groups",   1)
  inp.put_str("equation",        "sc")
  inp.put_int("inner_max_iters", 500)
  inp.put_dbl("inner_tolerance", 1e-4)
  inp.put_int("inner_print_out", 0)
  inp.put_int("outer_max_iters", 100)
  inp.put_dbl("outer_tolerance", 1e-4)
  inp.put_int("outer_print_out", 0)
  inp.put_str("bc_left",         "vacuum")
  inp.put_str("bc_right",        "reflect")
  inp.put_str("bc_bottom",       "reflect")
  inp.put_str("bc_top",          "vacuum")
  #-----------------------------------------------------------------------------#
  # MATERIAL
  #-----------------------------------------------------------------------------#
  mat = detran.Material.Create(1, 1, False)
  #
  mat.set_sigma_t(0, 0,    1.0)
  mat.set_sigma_s(0, 0, 0, 0.75)
  #
  mat.finalize()
  #-----------------------------------------------------------------------------#
  # GEOMETRY
  #-----------------------------------------------------------------------------#
  cm = [0.0, 10.0]
  fm = [50]
  mt = [0]
  mesh = detran.Mesh2D.Create(fm, fm, cm, cm, mt)

  # TEST LOOP
  for nq in range(0, nquad) :
    t = np.zeros(ntime)
    quad = detran.UniformEqual.Create(nq+1, 2)
    nang[nq] = quad.number_angles()
    for nt in range(0, ntime) :
      state = detran.State.Create(inp, mesh, quad)
      q_e = detran.ConstantSource.Create(mesh, quad, 1, 1.0)
      q_f = detran.FissionSource.Create(state, mesh, mat)
      bound = detran.Boundary2D.Create(inp, mesh, quad)
      solver = detran.GaussSeidel2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)
      start = time.time()
      solver.solve()
      t[nt] = (time.time() - start)
    times[nq] = np.average(t)
    stds[nq]  = np.std(t)

print times
print stds
print nang

