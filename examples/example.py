# Note: This tutorial assumes the Python front end has been built.

# The problem we'll solve is a one group, homogeneous slab, with
# a uniform, isotropic source, and vacuum boundaries.

from detran import *

# Next, create an input database, and fill it with some basic entries.

inp = InputDB.Create("optional_db_name")
inp.put_int("number_groups",   1)
inp.put_str("equation",        "dd")
inp.put_int("inner_max_iters", 100)
inp.put_dbl("inner_tolerance", 1e-4)
inp.put_int("inner_print_out", 0)
inp.put_str("bc_west",         "vacuum")
inp.put_str("bc_east",         "vacuum")


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

mat = Material.Create(1, 1, "optional_material_name")
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


# Finally, we can do some setup.  We first make the solver and
# perform initial setup.  This gets all the base stuff built,
# including the quadrature that is defined based on user input.

solver = Fixed1D(inp, mat, mesh)
solver.setup()

# We get the quadrature, which we put into a constant source.
# This is isotropic and uniform over all groups and space. The
# quadrature is needed in general since sources can be discrete.
# For diffusion problems, quad can be omitted for applicable
# source types (including ConstantSource and IsotropicSource).
# The source is set, and the solver gets a final setup.
quad = solver.quadrature()
#                    (number_groups, mesh, strength, quad)
q_e = ConstantSource.Create(1, mesh, 1.0, quad)
solver.set_source(q_e)
solver.set_solver()

# Now, solving is easy.
solver.solve()

# For fun, we can plot the flux.  Note, since Detran uses C++
# vectors for fluxes, we cast to Numpy arrays for plotting. There
# may be a less copy-intensive solution, but for the largest data
# sets, Python plotting would likely not be the best answer; rather
# the Silo interface is suggested for use with Visit or Paraview.
phi = np.asarray(solver.state().phi(0)) # "0" is for group 0
plot_mesh_function(mesh, phi)

