# Note: This tutorial assumes the Python front end has been built.

# The problem we'll solve is a one group, homogeneous slab, with
# a uniform, isotropic source, and vacuum boundaries.

# We first import the required modules:


import numpy as np
from detran import *


# Next, create an input database, and fill it with some basic entries.

inp = InputDB.Create()
inp.put_int("number_groups",   1)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 100)
inp.put_dbl("inner_tolerance", 1e-4)
inp.put_int("inner_print_out", 0)
inp.put_str("bc_left",         "vacuum")
inp.put_str("bc_right",        "vacuum")


# Now, we need the mesh.  The calling sequence gives the number
# of fine meshes in each coarse mesh, followed by the coarse mesh
# edges and the coarse mesh material assignment.  This sort of
# definition was inspired (for better or worse) from my days
# playing with DANTSYS.  In this example, we have just one
# coarse mesh for the slab, but you can see how this would
# simplify material and mesh assignment for less trivial problems.

mesh = Mesh1D.Create([100], [0.0, 10.0], [0])


# For the materials, we use take the total cross section to 
# be 1.0, and the scattering ratio to be 0.5.  The calling
# sequence is the number of groups, materials, and whether
# or not upscatter should be turned off explicitly (since
# the Material class defines group loop bounds).  The required
# call
# to `finalize` does things like adjust group bounds for
# outer iterations, but isn't too relevant here.

mat = Material.Create(1, 1, False)
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.5)
mat.finalize()


# Notice at this point the use of `Create` functions.  These are
# used because detran almost exclusively uses smart pointers
# for objects.  Since SWIG is used directly on the C++ source,
# a simple call like `Material(1, 1, False)` would produce
# a `Material` object, not a smart pointer to it.  The `Create`
# functions have the same signature as constructors but 
# return a smart pointer.

# For the quadrature, we'll use Gauss-Legendre (with 32 angles).

quad = GaussLegendre.Create(32)


# Finally, we can do some setup.  Eventually, this stuff might 
# go into a manager utility, but for now, there's a lot of
# explicit construction.  The `State` class contains all
# the fluxes, the eigenvalue (if applicable), and any other
# data that defines the solution.  The `ConstantSource`
# class is really a utility class that puts a constant
# source in all groups across all space.  Finally,
# the `Boundary1D` class and its 2D and 3D relatives
# contain the boundary fluxes.  In some sense, this
# could be viewed as an extension of the `State`.  The
# solver is a standard Gauss-Seidel solver with Source 
# Iteration for the inners.  Since this a one group 
# problem, the outer solver does essentially nothing.

state = State.Create(inp, mesh, quad)
q_e = ConstantSource.Create(mesh, quad, 1, 1.0)
q_f = FissionSource.Create(state, mesh, mat)
bound = Boundary1D.Create(inp, mesh, quad)
solver = GaussSeidel1D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)


# Now, solving is easy.

solver.solve()


# For fun, we can plot the flux. 
import matplotlib.pyplot as plt
phi = state.phi(0)
x   = np.linspace(0.0, 10.0, 101)
x   = 0.5*x[1] + x[0:100]
plt.plot(x, phi, linewidth = 2)
plt.grid(True)
plt.xlabel('x [cm]')
plt.ylabel('$\phi(x)$')
plt.title('A Monoenergetic Slab')
plt.savefig('tutorial.png')
plt.show()
