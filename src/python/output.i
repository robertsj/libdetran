
// File: index.xml

// File: classdetran_1_1Acceleration.xml
%feature("docstring") detran::Acceleration "

Base class for coarse mesh acceleration schemes.

All anticipated acceleration schemes have several shared features.
These include + Computing coarse mesh reaction rates + Computing
functions of coarse mesh boundary fluxes + Solving a low order
equation on the coarse mesh

C++ includes: Acceleration.hh ";

%feature("docstring")  detran::Acceleration::Acceleration "home
robertsj Research detran source src transport Acceleration cc
detran::Acceleration< D >::Acceleration(SP_mesh mesh, SP_material
material, SP_quadrature quadrature)

Constructor.

Parameters:
-----------

mesh:  Mesh smart pointer

material:  Material smart pointer

quadrature:  Quadrature smart pointer ";

%feature("docstring")  detran::Acceleration::~Acceleration "virtual
detran::Acceleration< D >::~Acceleration()

Virtual destructor. ";

%feature("docstring")  detran::Acceleration::initialize "virtual void
detran::Acceleration< D >::initialize(int level)=0

Create acceleration mesh given coarseness level and other setup.

By default, this initializes the coarse mesh by assigning a desired
number of fine meshes per coarse mesh. Extra meshes are assigned by
round-robin addition until all are assigned.

Clients may re-implement this to do more than just coarsen (e.g.
allocations).

Parameters:
-----------

level:  Desired number of fine meshes per coarse mesh ";

%feature("docstring")  detran::Acceleration::tally "virtual void
detran::Acceleration< D >::tally(int i, int j, int k, int o, int a,
face_flux_type psi)=0

Add contribution to an arbitrary function of the coarse mesh edge
flux.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

o:  octant

a:  angle within octant

psi:  edge angular flux ";

%feature("docstring")  detran::Acceleration::reset "virtual void
detran::Acceleration< D >::reset()=0

Reset for a new sweep. ";

%feature("docstring")  detran::Acceleration::set_group "void
detran::Acceleration< D >::set_group(int g) ";

%feature("docstring")  detran::Acceleration::fine_to_coarse "int
detran::Acceleration< D >::fine_to_coarse(int ijk, int dim) const

Homogenize the material data.

This function takes the current state vector and homogenizes the group
constants via flux-weighting.

Parameters:
-----------

state:  The current state vector

Get the coarse mesh index for a fine mesh

Parameters:
-----------

ijk:  fine mesh index

dim:  dimension of index

coarse mesh index ";

%feature("docstring")  detran::Acceleration::get_mesh "SP_mesh
detran::Acceleration< D >::get_mesh() const

Return the actual mesh. ";

%feature("docstring")  detran::Acceleration::get_coarse_mesh "SP_mesh
detran::Acceleration< D >::get_coarse_mesh() const

Return the coarse mesh. ";

%feature("docstring")  detran::Acceleration::get_material "SP_material detran::Acceleration< D >::get_material()

Return the actual material. ";

%feature("docstring")  detran::Acceleration::get_quadrature "SP_material detran::Acceleration< D >::get_quadrature()

Return the quadrature. ";

%feature("docstring")  detran::Acceleration::is_valid "bool
detran::Acceleration< D >::is_valid() const

Check if the object is in a valid state. ";


// File: classdetran__geometry_1_1Assembly.xml
%feature("docstring") detran_geometry::Assembly "

Simple square assembly.

Assembly represents a square array of pin cells that are assumed to
have the same meshing. The meshing is also assumed to be isotropic
(same in x and y) but not necessarily uniform.

C++ includes: Assembly.hh ";

%feature("docstring")  detran_geometry::Assembly::Assembly "home
robertsj Research detran source src geometry Assembly cc
detran_geometry::Assembly::Assembly(int dimension, vec_pincell
pincells, vec_int pincell_map)

Constructor.

Parameters:
-----------

pitch:  Pin cell pitch (assumed square)

radii:  Vector of fuel pin radii (can be zero length)

mat_map:  Region material map (cell-center outward)

meshes:  Number of evenly-spaced meshes per direction ";

%feature("docstring")  detran_geometry::Assembly::Assembly "detran_geometry::Assembly::Assembly(int dimension)

Constructor.

Parameters:
-----------

dimension:  Number of pins per row (e.g 17 for 17x17) ";

%feature("docstring")  detran_geometry::Assembly::mesh "Mesh::SP_mesh
detran_geometry::Assembly::mesh()

Return underlying meshed object. ";

%feature("docstring")  detran_geometry::Assembly::add_pincell "void
detran_geometry::Assembly::add_pincell(SP_pincell pin)

Add a pincell. ";

%feature("docstring")  detran_geometry::Assembly::finalize "void
detran_geometry::Assembly::finalize(vec_int pincell_map)

Mesh the assembly. ";

%feature("docstring")  detran_geometry::Assembly::dimension "int
detran_geometry::Assembly::dimension() const

Get dimension. ";

%feature("docstring")  detran_geometry::Assembly::number_pincells "int detran_geometry::Assembly::number_pincells() ";


// File: classdetran_1_1BoundaryBase.xml
%feature("docstring") detran::BoundaryBase "

Base class for boundary flux containers.

C++ includes: BoundaryBase.hh ";

%feature("docstring")  detran::BoundaryBase::BoundaryBase "detran::BoundaryBase< D >::BoundaryBase(SP_input input, SP_mesh mesh)

Constructor.

Parameters:
-----------

input:

mesh:  ";

%feature("docstring")  detran::BoundaryBase::~BoundaryBase "virtual
detran::BoundaryBase< D >::~BoundaryBase()

Virtual destructor. ";

%feature("docstring")  detran::BoundaryBase::set "virtual void
detran::BoundaryBase< D >::set(const size_t g)=0

Set the boundaries for a within-group solve.

This sets any boundaries that must be fixed for a solve, such as any
external boundary source.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryBase::update "virtual void
detran::BoundaryBase< D >::update(const size_t g)=0

Update the boundaries for each sweep.

This updates all incident boundary fluxes using the current outgoing
boundary fluxes in a group. What happens in the update is a function
of the underlying boundary condition.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryBase::update "virtual void
detran::BoundaryBase< D >::update(const size_t g, const size_t o,
const size_t a)=0

Update the boundaries for a single angle.

This is an alternative update that only updates the incident boundary
flux for a particular angle. When called within a sweep, this allows
the most recent boundary fluxes to be used, yielding a better
iteration.

This cannot be used for Krylov solvers.

Parameters:
-----------

g:  Group of current solve

o:  Octant being swept

a:  Angle within octant being swept ";

%feature("docstring")  detran::BoundaryBase::clear "virtual void
detran::BoundaryBase< D >::clear(const size_t g)=0

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryBase::psi "virtual void
detran::BoundaryBase< D >::psi(const size_t g, double *v, const int
inout, const int gs, bool onlyref=true)=0

Set the entire group boundary flux for reflecting sides.

This is in support of Krylov solvers.

Parameters:
-----------

g:  Group of current sweep

v:  Pointer to data used in Krylov solve ";

%feature("docstring")  detran::BoundaryBase::clear "void
detran::BoundaryBase< D >::clear()

Clear all groups. ";

%feature("docstring")  detran::BoundaryBase::has_reflective "bool
detran::BoundaryBase< D >::has_reflective() const

Does the boundary have any reflective conditions? ";

%feature("docstring")  detran::BoundaryBase::is_reflective "bool
detran::BoundaryBase< D >::is_reflective(const size_t side) const

Is a side reflective? ";

%feature("docstring")  detran::BoundaryBase::has_vacuum "bool
detran::BoundaryBase< D >::has_vacuum() const

Does the boundary have any vacuum conditions? ";

%feature("docstring")  detran::BoundaryBase::boundary_flux_size "size_t detran::BoundaryBase< D >::boundary_flux_size(const size_t
side) const

Number of boundary flux values in a group on a side. ";


// File: classdetran_1_1BoundaryCondition.xml
%feature("docstring") detran::BoundaryCondition "

Boundary condition for a surface.

C++ includes: BoundaryCondition.hh ";

%feature("docstring")  detran::BoundaryCondition::BoundaryCondition "detran::BoundaryCondition< D >::BoundaryCondition(Boundary_T
&boundary, const size_t side, SP_input input, SP_mesh mesh,
SP_quadrature quadrature) ";

%feature("docstring")  detran::BoundaryCondition::~BoundaryCondition "virtual detran::BoundaryCondition< D >::~BoundaryCondition()

Virtual destructor. ";

%feature("docstring")  detran::BoundaryCondition::set "virtual void
detran::BoundaryCondition< D >::set(const size_t g)=0

Set initial and/or fixed boundary condition. ";

%feature("docstring")  detran::BoundaryCondition::update "virtual
void detran::BoundaryCondition< D >::update(const size_t g)=0

Update a boundary following a sweep. ";

%feature("docstring")  detran::BoundaryCondition::update "virtual
void detran::BoundaryCondition< D >::update(const size_t g, const
size_t o, const size_t a)=0

Update a boundary for a given angle following a sweep. ";


// File: classdetran_1_1BoundaryConditionMOC.xml
%feature("docstring") detran::BoundaryConditionMOC "

Todo template on the Boundary type!!

C++ includes: BoundaryConditionMOC.hh ";

%feature("docstring")
detran::BoundaryConditionMOC::BoundaryConditionMOC "detran::BoundaryConditionMOC< D >::BoundaryConditionMOC(BoundaryMOC< D
> &boundary, const size_t side, SP_input input, SP_mesh mesh,
SP_quadrature quadrature) ";

%feature("docstring")
detran::BoundaryConditionMOC::~BoundaryConditionMOC "virtual
detran::BoundaryConditionMOC< D >::~BoundaryConditionMOC()

Virtual destructor. ";

%feature("docstring")  detran::BoundaryConditionMOC::set "virtual
void detran::BoundaryConditionMOC< D >::set(const size_t g)=0

Set initial and/or fixed boundary condition. ";

%feature("docstring")  detran::BoundaryConditionMOC::update "virtual
void detran::BoundaryConditionMOC< D >::update(const size_t g)=0

Update a boundary following a sweep. ";

%feature("docstring")  detran::BoundaryConditionMOC::update "virtual
void detran::BoundaryConditionMOC< D >::update(const size_t g, const
size_t o, const size_t a)=0

Update a boundary for a given angle following a sweep. ";


// File: classdetran_1_1BoundaryDiffusion.xml
%feature("docstring") detran::BoundaryDiffusion "

Container for diffusion boundary partial currents.

Unlike the boundary containers for transport, the diffusion boundary
is not responsible for boundary conditions. This is because these
conditions are built into the diffusion operator.

C++ includes: BoundaryDiffusion.hh ";

%feature("docstring")  detran::BoundaryDiffusion::BoundaryDiffusion "home robertsj Research detran source src boundary BoundaryDiffusion cc
detran::BoundaryDiffusion< D >::BoundaryDiffusion(SP_input input,
SP_mesh mesh)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh. ";

%feature("docstring")  detran::BoundaryDiffusion::set "void
detran::BoundaryDiffusion< D >::set(const size_t g)

Set the boundaries for a within-group solve.

This sets any boundaries that must be fixed for a solve, such as any
external boundary source.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryDiffusion::update "void
detran::BoundaryDiffusion< D >::update(const size_t g)

Update the boundaries for each sweep.

This updates all incident boundary fluxes using the current outgoing
boundary fluxes in a group. What happens in the update is a function
of the underlying boundary condition.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryDiffusion::update "void
detran::BoundaryDiffusion< D >::update(const size_t g, const size_t o,
const size_t a)

Update the boundaries for a single angle.

This is an alternative update that only updates the incident boundary
flux for a particular angle. When called within a sweep, this allows
the most recent boundary fluxes to be used, in effect producing a
Gauss-Seidel iteration.

This cannot be used for Krylov solvers.

Parameters:
-----------

g:  Group of current solve

o:  Octant being swept

a:  Angle within octant being swept ";

%feature("docstring")  detran::BoundaryDiffusion::clear "void
detran::BoundaryDiffusion< D >::clear(const size_t g)

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryDiffusion::psi "void
detran::BoundaryDiffusion< D >::psi(const size_t g, double *v, const
int inout, const int gs, bool onlyref=true)

Set the entire group boundary flux for reflecting sides.

This is is support of Krylov solvers.

Parameters:
-----------

g:  Group of current sweep

v:  Pointer to data used in Krylov solve ";

%feature("docstring")  detran::BoundaryDiffusion::get_input "SP_input
detran::BoundaryDiffusion< D >::get_input() const

Return the input. ";

%feature("docstring")  detran::BoundaryDiffusion::get_mesh "SP_mesh
detran::BoundaryDiffusion< D >::get_mesh() const

Return the mesh. ";

%feature("docstring")  detran::BoundaryDiffusion::clear "void
detran::BoundaryDiffusion< _3D >::clear(const size_t g)

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryDiffusion::clear "void
detran::BoundaryDiffusion< _2D >::clear(const size_t g)

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundaryDiffusion::clear "void
detran::BoundaryDiffusion< _1D >::clear(const size_t g)

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";


// File: classdetran_1_1BoundaryMOC.xml
%feature("docstring") detran::BoundaryMOC "

Boundary flux container for MOC problems.

The method of characteristics solves the transport equation by
sweeping along fixed tracks crossing the domain. Tracks begin and end
at a global boundary. This class stores the angular flux for each
track.

C++ includes: BoundaryMOC.hh ";

%feature("docstring")  detran::BoundaryMOC::BoundaryMOC "home
robertsj Research detran source src boundary BoundaryMOC cc home
robertsj Research detran source src boundary BoundaryMOC cc
detran::BoundaryMOC< D >::BoundaryMOC(SP_input input, SP_mesh mesh,
SP_quadrature quadrature)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

quadrature:  Angular quadrature. ";

%feature("docstring")  detran::BoundaryMOC::set "void
detran::BoundaryMOC< D >::set(const size_t g)

Set the boundaries for a within-group solve. ";

%feature("docstring")  detran::BoundaryMOC::update "void
detran::BoundaryMOC< D >::update(const size_t g)

Update the boundaries for each sweep. ";

%feature("docstring")  detran::BoundaryMOC::update "void
detran::BoundaryMOC< D >::update(const size_t g, const size_t o, const
size_t a)

brief Update the boundaries for a single angle. ";

%feature("docstring")  detran::BoundaryMOC::clear "void
detran::BoundaryMOC< D >::clear(const size_t g)

Clear the boundary container for a group. ";

%feature("docstring")  detran::BoundaryMOC::psi "void
detran::BoundaryMOC< D >::psi(const size_t g, double *v, const int
inout, const int gs, bool onlyref=true)

Set the entire group boundary flux for reflecting sides. ";

%feature("docstring")  detran::BoundaryMOC::feed_into "void
detran::BoundaryMOC< D >::feed_into(const size_t o1, const size_t a1,
const size_t t1, size_t &o2, size_t &a2, size_t &t2)

Assign the octant, angle, and track fed by an octant, angle, and
track.

To make the boundary condition as independent as possible, BoundaryMOC
knows what to do with the cardinal angle within an octant (i.e. with
polar implicit)

Parameters:
-----------

o1:  Incident octant

a1:  Incident angle within octant

t1:  Incident track

o2:  Outgoing octant

a2:  Outgoing angle within octant

t2:  Outgoing track ";

%feature("docstring")  detran::BoundaryMOC::feed_from "void
detran::BoundaryMOC< D >::feed_from(const size_t o1, const size_t a1,
const size_t t1, size_t &o2, size_t &a2, size_t &t2)

Assign the octant, angle, and track feeding an octant, angle, and
track.

Parameters:
-----------

o1:  Outgoing octant

a1:  Outgoing angle within octant

t1:  Outgoing track

o2:  Incident octant

a2:  Incident angle within octant

t2:  Incident track ";

%feature("docstring")  detran::BoundaryMOC::side_indices "const
vec2_int& detran::BoundaryMOC< D >::side_indices(const size_t side)
const

Return vector of octant, azimuth, track triplets for a side. ";

%feature("docstring")  detran::BoundaryMOC::BoundaryMOC "detran::BoundaryMOC< _1D >::BoundaryMOC(SP_input input, SP_mesh mesh,
SP_quadrature quadrature) ";

%feature("docstring")  detran::BoundaryMOC::BoundaryMOC "detran::BoundaryMOC< _3D >::BoundaryMOC(SP_input input, SP_mesh mesh,
SP_quadrature quadrature) ";


// File: classdetran_1_1BoundarySN.xml
%feature("docstring") detran::BoundarySN "

Boundary flux container for SN problems.

C++ includes: BoundarySN.hh ";

%feature("docstring")  detran::BoundarySN::BoundarySN "home robertsj
Research detran source src boundary BoundarySN cc home robertsj
Research detran source src boundary BoundarySN cc detran::BoundarySN<
D >::BoundarySN(SP_input input, SP_mesh mesh, SP_quadrature
quadrature)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

quadrature:  Angular quadrature. ";

%feature("docstring")  detran::BoundarySN::set "void
detran::BoundarySN< D >::set(const size_t g)

Set the boundaries for a within-group solve. ";

%feature("docstring")  detran::BoundarySN::update "void
detran::BoundarySN< D >::update(const size_t g)

Update the boundaries for each sweep. ";

%feature("docstring")  detran::BoundarySN::update "void
detran::BoundarySN< D >::update(const size_t g, const size_t o, const
size_t a)

Update the boundaries for a single angle. ";

%feature("docstring")  detran::BoundarySN::clear "void
detran::BoundarySN< D >::clear(const size_t g)

Clear the boundary container for a group. ";

%feature("docstring")  detran::BoundarySN::psi "void
detran::BoundarySN< D >::psi(const size_t g, double *v, const int
inout, const int gs, bool onlyref=true)

Set the entire group boundary flux for reflecting sides. ";

%feature("docstring")  detran::BoundarySN::get_input "SP_input
detran::BoundarySN< D >::get_input() const

Return the input. ";

%feature("docstring")  detran::BoundarySN::get_mesh "SP_mesh
detran::BoundarySN< D >::get_mesh() const

Return the mesh. ";

%feature("docstring")  detran::BoundarySN::get_quadrature "SP_quadrature detran::BoundarySN< D >::get_quadrature() const

Return the quadrature. ";

%feature("docstring")  detran::BoundarySN::clear "void
detran::BoundarySN< _2D >::clear(const size_t g)

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundarySN::clear "void
detran::BoundarySN< _1D >::clear(const size_t g)

Clear the boundary container for a group.

In some cases, a client might require homogeneous boundaries, perhaps
after a fixed boundary has been used to construct a right hand side
for a Krylov solve.

Parameters:
-----------

g:  Group of current solve ";

%feature("docstring")  detran::BoundarySN::psi "void
detran::BoundarySN< _3D >::psi(const size_t g, double *v, const int
inout, const int gs, bool onlyref)

Set the entire group boundary flux for reflecting sides.

This is in support of Krylov solvers.

Parameters:
-----------

g:  Group of current sweep

v:  Pointer to data used in Krylov solve ";

%feature("docstring")  detran::BoundarySN::psi "void
detran::BoundarySN< _2D >::psi(const size_t g, double *v, const int
inout, const int gs, bool onlyref)

Set the entire group boundary flux for reflecting sides.

This is in support of Krylov solvers.

Parameters:
-----------

g:  Group of current sweep

v:  Pointer to data used in Krylov solve ";

%feature("docstring")  detran::BoundarySN::psi "void
detran::BoundarySN< _1D >::psi(const size_t g, double *v, const int
inout, const int gs, bool onlyref)

Set the entire group boundary flux for reflecting sides.

This is in support of Krylov solvers.

Parameters:
-----------

g:  Group of current sweep

v:  Pointer to data used in Krylov solve ";


// File: classdetran__external__source_1_1BoundarySource.xml
%feature("docstring") detran_external_source::BoundarySource "

Base boundary source class.

C++ includes: BoundarySource.hh ";


// File: classdetran_1_1BoundaryTally.xml
%feature("docstring") detran::BoundaryTally "

Base class for recording angular flux function at coarse mesh
boundaries.

C++ includes: BoundaryTally.hh ";

%feature("docstring")  detran::BoundaryTally::BoundaryTally "detran::BoundaryTally< D >::BoundaryTally(SP_coarsemesh coarsemesh,
SP_quadrature quadrature, const size_t number_groups)

Constructor.

Parameters:
-----------

mesh:

quadrature:

number_groups:  ";

%feature("docstring")  detran::BoundaryTally::~BoundaryTally "virtual
detran::BoundaryTally< D >::~BoundaryTally()

Virtual destructor. ";

%feature("docstring")  detran::BoundaryTally::tally "virtual void
detran::BoundaryTally< D >::tally(const size_t i, const size_t j,
const size_t k, const size_t g, const size_t o, const size_t a, const
face_flux_type psi)=0

Add angular flux to the current tally.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

g:  group index

o:  octant

a:  angle within octant

psi:  edge angular flux ";

%feature("docstring")  detran::BoundaryTally::tally "virtual void
detran::BoundaryTally< D >::tally(const size_t i, const size_t j,
const size_t k, const size_t g, const size_t o, const size_t a, const
size_t d, const double psi)=0

Add angular flux to the tally for a single incident direction.

This is used to sweep over incident boundary containers so that the
incident conditions can be avoided for the sweep.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

g:  group index

o:  octant

a:  angle within octant

d:  axis index for the incident flux

psi:  edge angular flux ";

%feature("docstring")  detran::BoundaryTally::display "virtual void
detran::BoundaryTally< D >::display()=0

Print all the partial currents (for debugging) ";

%feature("docstring")  detran::BoundaryTally::reset "virtual void
detran::BoundaryTally< D >::reset(const size_t group)=0

Reset a group. ";


// File: structdetran_1_1BoundaryTraits.xml
%feature("docstring") detran::BoundaryTraits "

Boundary traits to simplify type access.

For fine mesh discretizations, the boundary data is stored for each
surface in the form of a value (1-D), 1-D vector (2-D), or 2-D vector
(3-D). This would be used for mesh-based discretizations, as in SN and
diffusion.

C++ includes: BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryTraits_3_01__1D_01_4.xml
%feature("docstring") detran::BoundaryTraits< _1D > " C++ includes:
BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryTraits_3_01__2D_01_4.xml
%feature("docstring") detran::BoundaryTraits< _2D > " C++ includes:
BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryTraits_3_01__3D_01_4.xml
%feature("docstring") detran::BoundaryTraits< _3D > " C++ includes:
BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryValue.xml
%feature("docstring") detran::BoundaryValue "

Boundary access.

Because we employ a templated boundary type, it can sometimes
complicate simple access within an otherwise general algorithm. This
accessor returns a boundary element given boundary spatial indices.

C++ includes: BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryValue_3_01__1D_01_4.xml
%feature("docstring") detran::BoundaryValue< _1D > " C++ includes:
BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryValue_3_01__2D_01_4.xml
%feature("docstring") detran::BoundaryValue< _2D > " C++ includes:
BoundaryTraits.hh ";


// File: structdetran_1_1BoundaryValue_3_01__3D_01_4.xml
%feature("docstring") detran::BoundaryValue< _3D > " C++ includes:
BoundaryTraits.hh ";


// File: classdetran__angle_1_1ChebyshevDPN.xml
%feature("docstring") detran_angle::ChebyshevDPN "

Chebyshev-Double PN quadrature.

The ChebyshevDPN quadrature uses Gauss-Chebyshev quadrature over one
azimuthal quadrant and DPN quadrature over the polar range.

Relevant database parameters: quad_number_polar_octant -- number polar
angles per octant

quad_number_azimuth_octant -- number azimuths per octant

C++ includes: ChebyshevDPN.hh ";

%feature("docstring")  detran_angle::ChebyshevDPN::ChebyshevDPN "home
robertsj Research detran source src angle ChebyshevDPN cc home
robertsj Research detran source src angle ChebyshevDPN cc
detran_angle::ChebyshevDPN::ChebyshevDPN(const size_t dim, const
size_t na, const size_t np)

Constructor.

Parameters:
-----------

dim:  Dimension (2 or 3)

na:  Number of azimuths per octant

np:  Number of polar angles per octant

name:   Quadrature name ";

%feature("docstring")  detran_angle::ChebyshevDPN::~ChebyshevDPN "virtual detran_angle::ChebyshevDPN::~ChebyshevDPN()

Virtual destructor. ";


// File: classdetran__angle_1_1ChebyshevLegendre.xml
%feature("docstring") detran_angle::ChebyshevLegendre "

Chebyshev-Legendre quadrature.

The ChebyshevLegendre quadrature uses Gauss-Chebyshev quadrature over
one azimuthal quadrant and Gauss-Legendre quadrature over the polar
range.

Relevant database parameters: quad_number_polar_octant -- number polar
angles per octant

quad_number_azimuth_octant -- number azimuths per octant

C++ includes: ChebyshevLegendre.hh ";

%feature("docstring")
detran_angle::ChebyshevLegendre::ChebyshevLegendre "home robertsj
Research detran source src angle ChebyshevLegendre cc
detran_angle::ChebyshevLegendre::ChebyshevLegendre(const size_t dim,
const size_t na, const size_t np)

Constructor.

Parameters:
-----------

dim:  Dimension (2 or 3)

na:  Number of azimuths per octant

np:  Number of polar angles per octant

name:   Quadrature name ";

%feature("docstring")
detran_angle::ChebyshevLegendre::~ChebyshevLegendre "virtual
detran_angle::ChebyshevLegendre::~ChebyshevLegendre()

Virtual destructor. ";


// File: classdetran_1_1CMR.xml
%feature("docstring") detran::CMR "

Coarse Mesh Rebalance.

Consider the one group problem \\\\[ \\\\Big ( \\\\mu
\\\\frac{\\\\partial}{\\\\partial x} + \\\\eta
\\\\frac{\\\\partial}{\\\\partial y} + \\\\xi
\\\\frac{\\\\partial}{\\\\partial z} \\\\Big ) \\\\psi +
\\\\Sigma_t(\\\\vec{r}) \\\\psi(\\\\vec{r}, \\\\hat{\\\\Omega}) =
\\\\Sigma_s(\\\\vec{r}) \\\\phi(\\\\vec{r}) + q(\\\\vec{r},
\\\\hat{\\\\Omega}) \\\\, . \\\\] where we have explicitly separated
the within-group scatter from the external source.

Suppose we integrate this over a coarse cell of volume $ V =
\\\\Delta_x \\\\Delta_y \\\\Delta_z $ and over the angular space, $
4\\\\pi $. The streaming terms give rise to terms like

\\\\[ \\\\begin{split} & \\\\int_{4\\\\pi} \\\\mu
\\\\int^{\\\\Delta_y}_{0} \\\\int^{\\\\Delta_z}_{0} dy dz \\\\Bigg (
\\\\frac{\\\\partial \\\\psi}{\\\\partial x} \\\\Big |_{x=\\\\Delta_x}
-\\\\frac{\\\\partial \\\\psi}{\\\\partial x} \\\\Big |_{x=0} \\\\Bigg
) = \\\\int^{\\\\Delta_y}_{0} \\\\int^{\\\\Delta_z}_{0} dy dz \\\\Big
( J_x(\\\\Delta_x, y, z) - J_x(0, y, z) \\\\Big ) \\\\, ,
\\\\end{split} \\\\]

where $ J_x $ is the $x$-directed partial current. Similar terms can
be found for the other directions.

We can write the integrated transport equation as

\\\\[ \\\\Delta_y \\\\Delta_z \\\\Big ( \\\\bar{J}_x(\\\\Delta_x) -
\\\\bar{J}_x(0) \\\\Big ) + \\\\Delta_x \\\\Delta_z \\\\Big (
\\\\bar{J}_y(\\\\Delta_y) - \\\\bar{J}_y(0) \\\\Big ) + \\\\Delta_x
\\\\Delta_y \\\\Big ( \\\\bar{J}_z(\\\\Delta_z) - \\\\bar{J}_z(0)
\\\\Big ) + \\\\Delta_x \\\\Delta_y \\\\Delta_z \\\\bar{\\\\Sigma_r}
\\\\bar{\\\\phi}(\\\\vec{r}) = \\\\Delta_x \\\\Delta_y \\\\Delta_z
\\\\bar{q}(\\\\vec{r}) \\\\]

where we have defined the average current across a face of the box as

\\\\[ \\\\bar{J}_x(x) = \\\\frac{\\\\int^{\\\\Delta_y}_{0}
\\\\int^{\\\\Delta_z}_{0} dy dz J_x(x, y, z)}
{\\\\int^{\\\\Delta_y}_{0} \\\\int^{\\\\Delta_z}_{0} dy dz} \\\\, .
\\\\]

and similar averages for the removal rate and source. Note, it is the
removal rate of interest on the left. This requires a coarse mesh
averaged removal cross section and scalar flux. The removal cross
section is simply the total cross section minus the within-group
scattering cross section. As a consequence, the source $ q $ does not
include the within-group scatter source.

We simplify notation by writing

\\\\[ \\\\sum_{p \\\\in (x,y,z)} \\\\frac{1}{\\\\Delta_{p}} \\\\Big (
\\\\bar{J}_{p}(\\\\Delta_{p}) - \\\\bar{J}_{p}(0) \\\\Big ) +
\\\\bar{\\\\Sigma_r} \\\\bar{\\\\phi}(\\\\vec{r}) =
\\\\bar{q}(\\\\vec{r}) \\\\]

This neutron balance equation is not satisfied when the flux in not
converged. However, if we can somehow force the flux to satisfy
balance at the coarse level, then we can apply any correction to the
fine level. Let us define rebalance factors $ f_i $ such that

\\\\[ \\\\sum_{p \\\\in (x,y,z)} \\\\frac{1}{\\\\Delta_{p}} \\\\Big (
\\\\bar{J}_{p}(\\\\Delta_{p}) - \\\\bar{J}_{p}(0) \\\\Big ) + f_i
\\\\bar{\\\\Sigma_r} \\\\bar{\\\\phi}_i = \\\\bar{q}(\\\\vec{r}) \\\\]

C++ includes: CMR.hh ";

%feature("docstring")  detran::CMR::CMR "home robertsj Research
detran source src transport CMR cc detran::CMR< D >::CMR(SP_input
input, SP_material material, SP_coarsemesh coarsemesh, SP_currenttally
currenttally)

Constructor.

Parameters:
-----------

input:  Input database

material:  Material database

coarsemesh:  Coarse mesh

currenttally:  Current tally ";


// File: classdetran_1_1CoarseMesh.xml
%feature("docstring") detran::CoarseMesh "C++ includes: CoarseMesh.hh
";

%feature("docstring")  detran::CoarseMesh::CoarseMesh "CoarseMesh::CoarseMesh(SP_mesh fine_mesh, const size_t level) ";

%feature("docstring")  detran::CoarseMesh::get_fine_mesh "SP_mesh
detran::CoarseMesh::get_fine_mesh() const

Get the original fine mesh. ";

%feature("docstring")  detran::CoarseMesh::get_coarse_mesh "SP_mesh
detran::CoarseMesh::get_coarse_mesh() const

Get the coarse mesh. ";

%feature("docstring")  detran::CoarseMesh::fine_to_coarse "size_t
detran::CoarseMesh::fine_to_coarse(const size_t ijk, const size_t dim)
const

Get the coarse mesh index for a fine mesh.

Parameters:
-----------

ijk:  fine mesh index

dim:  dimension of index

coarse mesh index ";

%feature("docstring")  detran::CoarseMesh::coarse_edge_flag "size_t
detran::CoarseMesh::coarse_edge_flag(const size_t ijk, const size_t
dim) const

Determine whether a fine mesh edge is on a coarse mesh edge.

Parameters:
-----------

ijk:  fine mesh edge index

dim:  dimension of index

nonnegative index of coarse edge (otherwise -1) ";


// File: classCoarseMesh.xml
%feature("docstring") CoarseMesh "

Create a coarse mesh for acceleration.

This is used for nonlinear acceleration schemes, and, in principle,
could be used to develop multigrid preconditioners for Krylov solvers.

This is a very limited interface. A possible second approach (for the
user) would be to assign a coarse mesh map (independent of the coarse
mesh used for material assignment)

C++ includes: CoarseMesh.hh ";


// File: classdetran__angle_1_1Collocated.xml
%feature("docstring") detran_angle::Collocated "

A MOC quadrature with a finite set of origins along a side.

This quadrature has a number of unique features. It satisfies cyclic
tracking requirements in a rectangular region. It also uses a finite
set of points along a surface as track entrance (and exit) points. In
this first implementation, we consider only square regions.

The user sets the number of spatial points along a side of a unit
cell. Currently, this must be a power of 3, and only 9 and 27 are
supported. Additionally, a multiple may be assigned so that the base
points are repeated on adjacent unit cells. This lets us define the
spacing in terms of a pin cell, which is then repeated on the assembly
level.

The angles (between $ \\\\pi/4 $ and $ \\\\pi/2 $) are uniquely
defined by \\\\[ \\\\tan{\\\\phi_i} = 3^i \\\\, , \\\\,\\\\,\\\\, i =
0, \\\\, \\\\cdots n \\\\, , \\\\] where \\\\[ n =
\\\\mathrm{floor}(\\\\log_3(N)+1) \\\\, \\\\] and N is the number of
spatial points. Note that n/N is maximized when N is a power of three.
Using 9 spatial points surface yields 5 angles per quadrant. 27 points
yields 7 angles per quadrant. The other angles are defined by symmetry
about $ \\\\pi/4 $.

For weights, the user can use an arc length approximation or a
quadrature rule defined for the points that may provide better
accuracy.

C++ includes: Collocated.hh ";

%feature("docstring")  detran_angle::Collocated::Collocated "home
robertsj Research detran source src angle Collocated cc
detran_angle::Collocated::Collocated(size_t dim, size_t
num_azimuths_octant, size_t multiplier, size_t num_polar, std::string
polar)

Constructor.

Parameters:
-----------

dim:  Problem dimension (only 2 supported for now)

num_azimuths_octant:  Number of azimuths per octant

num_space:  Number of tracks per azimuth

num_polar:  Number of polar angles in half space

polar:  Polar quadrature string identifier ";

%feature("docstring")  detran_angle::Collocated::~Collocated "detran_angle::Collocated::~Collocated() ";


// File: structdetran__ioutils_1_1compound__type.xml
%feature("docstring") detran_ioutils::compound_type "

Store contents of map key/value pair.

C++ includes: IO_HDF5_Traits.hh ";


// File: structdetran__ioutils_1_1compound__type__traits.xml
%feature("docstring") detran_ioutils::compound_type_traits "

Define the value type for the value of a map key/value pair.

To map std::map to HDF5, we employ a compound type that stores both
the key and value together. These traits help reduce code volume.

C++ includes: IO_HDF5_Traits.hh ";


// File: structdetran__ioutils_1_1compound__type__traits_3_01std_1_1string_01_4.xml
%feature("docstring") detran_ioutils::compound_type_traits<
std::string > " C++ includes: IO_HDF5_Traits.hh ";


// File: structdetran__ioutils_1_1compound__type__traits_3_01std_1_1vector_3_01double_01_4_01_4.xml
%feature("docstring") detran_ioutils::compound_type_traits<
std::vector< double > > " C++ includes: IO_HDF5_Traits.hh ";


// File: structdetran__ioutils_1_1compound__type__traits_3_01std_1_1vector_3_01int_01_4_01_4.xml
%feature("docstring") detran_ioutils::compound_type_traits<
std::vector< int > > " C++ includes: IO_HDF5_Traits.hh ";


// File: classdetran__external__source_1_1ConstantSource.xml
%feature("docstring") detran_external_source::ConstantSource "

Defines a single isotropic source everywhere.

C++ includes: ConstantSource.hh ";

%feature("docstring")
detran_external_source::ConstantSource::ConstantSource "home robertsj
Research detran source src external_source ConstantSource cc
detran_external_source::ConstantSource::ConstantSource(size_t
number_groups, SP_mesh mesh, double strength, SP_quadrature
quadrature=SP_quadrature(0))

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

strength:  Source strength in all groups and all space; the unit is n
/cc-sec

quadrature:  Pointer to quadrature (optional) ";

%feature("docstring")  detran_external_source::ConstantSource::source
"double detran_external_source::ConstantSource::source(const size_t
cell, const size_t group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran_external_source::ConstantSource::source
"double detran_external_source::ConstantSource::source(const size_t
cell, const size_t group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";


// File: classdetran__geometry_1_1Core.xml
%feature("docstring") detran_geometry::Core "

Simple 2-D core of assemblies.

C++ includes: Core.hh ";

%feature("docstring")  detran_geometry::Core::Core "home robertsj
Research detran source src geometry Core cc
detran_geometry::Core::Core(int dimension, vec_assembly assemblies,
vec_int assembly_map)

Constructor.

Parameters:
-----------

pitch:  Pin cell pitch (assumed square)

radii:  Vector of fuel pin radii (can be zero length)

mat_map:  Region material map (cell-center outward)

meshes:  Number of evenly-spaced meshes per direction ";

%feature("docstring")  detran_geometry::Core::Core "detran_geometry::Core::Core(int dimension)

Constructor.

Parameters:
-----------

dimension:  Number of pins per row (e.g 17 for 17x17) ";

%feature("docstring")  detran_geometry::Core::mesh "Mesh::SP_mesh
detran_geometry::Core::mesh()

Return underlying meshed object. ";

%feature("docstring")  detran_geometry::Core::add_assembly "void
detran_geometry::Core::add_assembly(SP_assembly assembly)

Add an assembly. ";

%feature("docstring")  detran_geometry::Core::finalize "void
detran_geometry::Core::finalize(vec_int assembly_map)

Mesh the assembly. ";

%feature("docstring")  detran_geometry::Core::pincell_index "int
detran_geometry::Core::pincell_index(int i, int j)

Pincell index. ";


// File: classdetran_1_1CurrentTally.xml
%feature("docstring") detran::CurrentTally "

Records partial currents through coarse mesh surfaces.

To illustrate, consider the leakage term of the transport equation in
one coarse cell:

\\\\[ \\\\Big ( \\\\mu \\\\frac{\\\\partial}{\\\\partial x} + \\\\eta
\\\\frac{\\\\partial}{\\\\partial y} + \\\\xi
\\\\frac{\\\\partial}{\\\\partial z} \\\\Big ) \\\\psi(x, \\\\mu,
\\\\eta, \\\\xi) \\\\, . \\\\]

Assuming a cell of volume $ V = \\\\Delta_x \\\\Delta_y \\\\Delta_z $,
we integrate the first term over the cell to get

\\\\[ \\\\begin{split} & \\\\mu \\\\int^{\\\\Delta_y}_{0}
\\\\int^{\\\\Delta_z}_{0} dy dz \\\\Bigg ( \\\\frac{\\\\partial
\\\\psi}{\\\\partial x} \\\\Big |_{x=\\\\Delta_x}
-\\\\frac{\\\\partial \\\\psi}{\\\\partial x} \\\\Big |_{x=0} \\\\Bigg
) = \\\\int^{\\\\Delta_y}_{0} \\\\int^{\\\\Delta_z}_{0} dy dz \\\\Big
( J_x(\\\\Delta_x, y, z) - J_x(0, y, z) \\\\Big ) \\\\, ,
\\\\end{split} \\\\]

where $ J_x $ is the $x$-directed partial current. Similar terms can
be found for the other directions.

In 2D and 3D cases, the current must be integrated in space, which is
done using a simple mid-point rule consistent with the underlying
spatial discretizations currently available in Detran. If higher order
methods are implemented, equation-dependent currents would be
required.

C++ includes: CurrentTally.hh ";

%feature("docstring")  detran::CurrentTally::CurrentTally "detran::CurrentTally< D >::CurrentTally(SP_coarsemesh mesh,
SP_quadrature quadrature, const size_t number_groups)

Constructor.

Parameters:
-----------

mesh:

quadrature:

number_groups:  ";

%feature("docstring")  detran::CurrentTally::tally "void
detran::CurrentTally< D >::tally(const size_t i, const size_t j, const
size_t k, const size_t g, const size_t o, const size_t a, const
face_flux_type psi)

Add angular flux to the current tally.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

g:  group index

o:  octant

a:  angle within octant

psi:  edge angular flux ";

%feature("docstring")  detran::CurrentTally::tally "void
detran::CurrentTally< D >::tally(const size_t i, const size_t j, const
size_t k, const size_t g, const size_t o, const size_t a, const size_t
d, const double psi)

Add angular flux to the tally for a single incident direction.

This is used to sweep over incident boundary containers so that the
incident conditions can be avoided for the sweep.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

g:  group index

o:  octant

a:  angle within octant

d:  axis index for the incident flux

psi:  edge angular flux ";

%feature("docstring")  detran::CurrentTally::partial_current "double
detran::CurrentTally< D >::partial_current(const size_t i, const
size_t j, const size_t k, const size_t g, const size_t axis, const
size_t sense)

Get the partial current from a surface and sense.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

g:  group index

axis:  0, 1, or 2 for x, y, or z

sense:  true for positive (e.g. +x) ";

%feature("docstring")  detran::CurrentTally::display "void
detran::CurrentTally< D >::display()

Print all the partial currents (for debugging) ";

%feature("docstring")  detran::CurrentTally::reset "void
detran::CurrentTally< D >::reset(const size_t group)

Reset a group. ";

%feature("docstring")  detran::CurrentTally::tally "void
detran::CurrentTally< _1D >::tally(const size_t i, const size_t j,
const size_t k, const size_t g, const size_t o, const size_t a, const
face_flux_type psi)

Add angular flux to the current tally.

Parameters:
-----------

i:  x mesh index

j:  y mesh index

k:  z mesh index

g:  group index

o:  octant

a:  angle within octant

psi:  edge angular flux ";


// File: classdetran_1_1DiffusionEigensolver.xml
%feature("docstring") detran::DiffusionEigensolver "

Eigensolver for multigroup diffusion equations.

C++ includes: DiffusionEigensolver.hh ";

%feature("docstring")
detran::DiffusionEigensolver::DiffusionEigensolver "home robertsj
Research detran source src solvers eigen DiffusionEigensolver cc
detran::DiffusionEigensolver< D >::DiffusionEigensolver(SP_input
input, SP_material material, SP_mesh mesh, SP_state state)

Constructor.

Parameters:
-----------

input:  parameter database

material:  material database

mesh:  mesh definition

state:  state vector ";

%feature("docstring")  detran::DiffusionEigensolver::solve "void
detran::DiffusionEigensolver< D >::solve()

Solve the system. ";


// File: classdetran_1_1DiffusionGainOperator.xml
%feature("docstring") detran::DiffusionGainOperator "C++ includes:
DiffusionGainOperator.hh ";

%feature("docstring")
detran::DiffusionGainOperator::DiffusionGainOperator "home robertsj
Research detran source src solvers mg DiffusionGainOperator cc
detran::DiffusionGainOperator::DiffusionGainOperator(SP_input input,
SP_material material, SP_mesh mesh, bool adjoint=false)

Constructor.

Parameters:
-----------

input:  parameter database

material:  material database

mesh:  mesh definition ";


// File: classdetran_1_1DiffusionLossOperator.xml
%feature("docstring") detran::DiffusionLossOperator "

Loss operator for multigroup diffusion problems.

The loss operator represents neutron losses through interactions and
leakage from a volume. In Detran, we employ a mesh-centered finite
difference approximation. Moreover, we use an energy block structure
so that the diffusion operators match with the underlying transport
structure, which is useful for preconditioning applications.

The loss operator can be constructed with or without a fission
contribution. For multiplying fixed source problems, it is numerically
more efficient to account for fission implicitly by bringing it to the
left hand side. However, there are cases where performing fission
iteration is warranted, in which case not including the fission source
is required.

C++ includes: DiffusionLossOperator.hh ";

%feature("docstring")
detran::DiffusionLossOperator::DiffusionLossOperator "home robertsj
Research detran source src solvers mg DiffusionLossOperator cc
detran::DiffusionLossOperator::DiffusionLossOperator(SP_input input,
SP_material material, SP_mesh mesh, const bool include_fission, const
size_t cutoff=0, const bool adjoint=false, const double keff=1.0)

Constructor.

Parameters:
-----------

input:  Pointer to input parameters

material:  Pointer to materials

mesh:  Pointer to mesh

include_fission:  Flag for including fission source implicitly

cutoff:  Lowest group to include in the operator

adjoint:  Adjoint flag

keff:  Fission scaling factor ";

%feature("docstring")  detran::DiffusionLossOperator::construct "void
detran::DiffusionLossOperator::construct(double keff=1.0)

Rebuild the matrix for a new fission scaling constant.

This allows the client to rebuild the matrix after initial
construction. This is useful for response function generation as a
function of keff or for time-dependent problems in which the pseudo-
coefficients changes with time.

Parameters:
-----------

keff:  Scaling parameter for fission source ";


// File: classdetran_1_1DimensionTraits.xml
%feature("docstring") detran::DimensionTraits "

Allows templates using dimension without <int n>=\"\"> format.

There may eventually be further uses for the class.

C++ includes: DimensionTraits.hh ";


// File: classdetran__external__source_1_1DiscreteSource.xml
%feature("docstring") detran_external_source::DiscreteSource "

Discrete source definition.

C++ includes: DiscreteSource.hh ";

%feature("docstring")
detran_external_source::DiscreteSource::DiscreteSource "home robertsj
Research detran source src external_source DiscreteSource cc
detran_external_source::DiscreteSource::DiscreteSource(size_t
number_groups, SP_mesh mesh, vec3_dbl spectra, vec_int map,
SP_quadrature quadrature)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

spectra:  Angle-dependent spectra, size = [#spectra][#groups][#angle]

map:  Map of where spectra are located, size = [#cells]

quadrature:  Quadrature ";

%feature("docstring")  detran_external_source::DiscreteSource::source
"double detran_external_source::DiscreteSource::source(const size_t
cell, const size_t group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran_external_source::DiscreteSource::source
"double detran_external_source::DiscreteSource::source(const size_t
cell, const size_t group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";


// File: classdetran__angle_1_1DPN.xml
%feature("docstring") detran_angle::DPN "

Implements the double GaussLegendre quadrature.

This is identical to Gauss-Legendre quadrature except that the
quadrature parameters are based on the integral \\\\[ \\\\int^{1}_{0}
f(x) dx \\\\approx \\\\sum^m_{i=1} W_m f(x_i) \\\\, . \\\\] Doing this
effectively treats the left and right half spaces independently, which
can work well for functions that have completely different behavior in
each half space but within those half spaces, the behavior is
polynomial-like.

The resulting abscissa are those of the GL quadrature but scaled and
shifted.

Relevant database parameters: quad_number_polar_octant -- number of
abscissa per half space

C++ includes: DPN.hh ";

%feature("docstring")  detran_angle::DPN::DPN "home robertsj Research
detran source src angle DPN cc detran_angle::DPN::DPN(const size_t
number_polar_octant)

Constructor.

Parameters:
-----------

number_polar_octant:  Number of polar angles per octant ";


// File: classdetran__angle_1_1DTN.xml
%feature("docstring") detran_angle::DTN "

Implements the double GaussChebyshev quadrature.

This is identical to Gauss-Chebyshev quadrature except that the
quadrature parameters are based on the integral \\\\[ \\\\int^{1}_{0}
f(x) dx \\\\approx \\\\sum^m_{i=1} W_m f(x_i) \\\\, . \\\\] Doing this
effectively treats the left and right half spaces independently, which
can work well for functions that have completely different behavior in
each half space but within those half spaces, the behavior is
polynomial-like.

The resulting abscissa are those of the GC quadrature but scaled and
shifted.

Relevant database parameters: quad_number_polar_octant -- number of
abscissa per half space

gausschebyshev_normalize -- normalize half space weight to 1 [false]

C++ includes: DTN.hh ";

%feature("docstring")  detran_angle::DTN::DTN "home robertsj Research
detran source src angle DTN cc home robertsj Research detran source
src angle DTN cc detran_angle::DTN::DTN(const size_t
number_polar_octant, const bool normalize=false)

Constructor.

Parameters:
-----------

number_polar_octant:  number of angles per half space

normalize:  normalize half space weights to 1 ";


// File: classdetran_1_1EigenArnoldi.xml
%feature("docstring") detran::EigenArnoldi "

Solves the eigenvalue problem using the Arnoldi method.

The eigenvalue problem can be cast in the form \\\\[ \\\\mathbf{A}d =
kd \\\\, \\\\] where $ d $ is the fission density and $ k $ is the
eigenvalue. See Eigensolver for more details on this formulation.

As mentioned in Eigensolver, Krylov methods require an adequately-
converged multigroup solve. The paper by Warsa et al. talks more about
this, and it seems the multigroup convergence is most important in the
first several iterations. Thus, it might be worth implementing a
dynamic tolerance at some point.

Detran must be configured with SLEPc for Arnoldi to be available via
the callow interface. Otherwise, the callow implementation of power
iteration is used.

C++ includes: EigenArnoldi.hh ";

%feature("docstring")  detran::EigenArnoldi::EigenArnoldi "home
robertsj Research detran source src solvers eigen EigenArnoldi cc
detran::EigenArnoldi< D >::EigenArnoldi(SP_mg_solver mg_solver)

Constructor.

Parameters:
-----------

mg_solver:  Multigroup solver ";

%feature("docstring")  detran::EigenArnoldi::solve "void
detran::EigenArnoldi< D >::solve()

Solve the eigenvalue problem. ";


// File: classdetran_1_1EigenDiffusion.xml
%feature("docstring") detran::EigenDiffusion "

Solves the diffusion eigenvalue directly.

While the other eigensolvers can be used in conjuction with the
multigroup diffusion solver, in this eigensolver, both the loss and
gains operator are explicitly constructed for solution via a
generalized eigensolver. Including this capability is mostly for
benchmarking the performance of the other approach.

C++ includes: EigenDiffusion.hh ";

%feature("docstring")  detran::EigenDiffusion::EigenDiffusion "home
robertsj Research detran source src solvers eigen EigenDiffusion cc
detran::EigenDiffusion< D >::EigenDiffusion(SP_mg_solver mg_solver)

Constructor.

Parameters:
-----------

mg_solver:  Multigroup solver ";

%feature("docstring")  detran::EigenDiffusion::solve "void
detran::EigenDiffusion< D >::solve()

Solve the eigenvalue problem. ";


// File: classdetran_1_1EigenPI.xml
%feature("docstring") detran::EigenPI "

Solves the eigenvalue problem via the power method.

The eigenvalue problem can be cast in the form \\\\[ \\\\mathbf{A}d =
kd \\\\, \\\\] where $ d $ is the fission density and $ k $ is the
eigenvalue. See Eigensolver for more details on this formulation.

The power method solves the eigenproblem using the iteration \\\\[
d^{l+1} \\\\leftarrow \\\\mathbf{A} d^{l} / k^{l} \\\\] and \\\\[
k^{l+1} \\\\leftarrow || d^{l+1} || \\\\, . \\\\] Traditionally, the
norm used to define the updated eigenvalue is a fission-weighted sum
of the group fluxes, which is equivalent to an L1 norm of the density
if all the fluxes and fission cross sections are positive, true for
all physical problems (and barring numerical issues due to
discretization). Here, we use the L1 norm, and begin with $ || d^{0}
|| = 1 $.

Note, this is a hand-coded power iteration implementation that can be
used with nonlinear acceleration.

C++ includes: EigenPI.hh ";

%feature("docstring")  detran::EigenPI::EigenPI "home robertsj
Research detran source src solvers eigen EigenPI cc detran::EigenPI< D
>::EigenPI(SP_mg_solver mg_solver)

Constructor.

Parameters:
-----------

mg_solver:  Multigroup solver ";

%feature("docstring")  detran::EigenPI::solve "void detran::EigenPI<
D >::solve()

Solve the eigenvalue problem. ";


// File: classcallow_1_1EigenSolver.xml
%feature("docstring") callow::EigenSolver "

Base class for iterative eigensolvers.

We solve generalized eigenvalue problems of the form \\\\[
\\\\mathbf{A}x = \\\\lambda \\\\mathbf{B} x \\\\] via iterative
methods. For the case when $ \\\\mathbf{B} \\\\ne \\\\mathbf{I} $ we
need to solve \\\\[ \\\\mathbf{B}^{-1}\\\\mathbf{A}x = \\\\lambda x
\\\\] A system is \"solved\" when the norm of the residual is small
enough or some maximum iteration count is reached. The residual is
defined \\\\[ \\\\mathbf{A}x - \\\\lambda \\\\mathbf{B} x \\\\] and
its norm is \\\\[ r = || (\\\\mathbf{A} - \\\\lambda \\\\mathbf{B})x||
\\\\] By default, the L2 norm is used, though L1 and Linf are also
recorded can can be used. In some cases, a different norm is
warranted, perhaps based on physics. This can be implemented by
derived classes.

Currently, we implement only the power method, inverse iteration, and
Rayleigh quotient. However, other solvers are available if SLEPc is
enabled. Additionally, our structure is really intended for the
dominant mode. Other modes will require different handling.

C++ includes: EigenSolver.hh ";

%feature("docstring")  callow::EigenSolver::EigenSolver "home
robertsj Research detran source src callow solver EigenSolver cc
callow::EigenSolver::EigenSolver(const double tol=1e-6, const int
maxit=100, std::string name=\"solver\") ";

%feature("docstring")  callow::EigenSolver::~EigenSolver "virtual
callow::EigenSolver::~EigenSolver() ";

%feature("docstring")  callow::EigenSolver::set_operators "void
callow::EigenSolver::set_operators(SP_matrix A, SP_matrix
B=SP_matrix(0), SP_db db=SP_db(0))

Sets the operators for the problem.

This allows for the system

Parameters:
-----------

A:  left side operator

B:  optional right side operator (to be inverted)

db:  optional database for solver and preconditioner options ";

%feature("docstring")  callow::EigenSolver::set_tolerances "void
callow::EigenSolver::set_tolerances(const double tol, const int maxit)

Set the convergence criteria.

Parameters:
-----------

tol:  m tolerance (||r_n|| < tol)

maxit:  maximum iterations (n < maxit) ";

%feature("docstring")  callow::EigenSolver::set_monitor_level "void
callow::EigenSolver::set_monitor_level(const int v)

Set the diagnostic level.

Higher levels means more output.

Parameters:
-----------

v:  diagnostic level ";

%feature("docstring")  callow::EigenSolver::solve "int
callow::EigenSolver::solve(Vector &x, Vector &x0)

Solve the eigenvalue problem.

Upon return, the initial guess is *not* guaranteed to be unchanged, as
it may be used as temporary storage.

Parameters:
-----------

x:  vector to fill with solution

x0:  initial guess ";

%feature("docstring")  callow::EigenSolver::solve "int
callow::EigenSolver::solve(SP_vector x, SP_vector x0)

Solve the eigenvalue problem (SP variant) ";

%feature("docstring")  callow::EigenSolver::residual_norms "std::vector<double> callow::EigenSolver::residual_norms()

Return the residual norms. ";

%feature("docstring")  callow::EigenSolver::A "SP_matrix
callow::EigenSolver::A()

Get the left operator. ";

%feature("docstring")  callow::EigenSolver::B "SP_matrix
callow::EigenSolver::B()

Get the right operator. ";

%feature("docstring")  callow::EigenSolver::linearsolver "SP_linearsolver callow::EigenSolver::linearsolver()

Get the linear solver. ";

%feature("docstring")  callow::EigenSolver::eigenvalue "double
callow::EigenSolver::eigenvalue()

Get the eigenvalue. ";


// File: classdetran_1_1Eigensolver.xml
%feature("docstring") detran::Eigensolver "

Base class for solving the eigenvalue problem.

The Eigenvalue Problem The steady-state balance of neutrons in a
fissile system is characterized by the multigroup eigenvalue problem

\\\\[ @hat{\\\\Omega} \\\\cdot \\\\nabla @psi +
@Sigma_t(\\\\vec{r},E_g) \\\\psi(\\\\vec{r},\\\\hat{\\\\Omega},E_g) =
@frac{1}{4\\\\pi} \\\\sum^G_{g'} \\\\int_{4\\\\pi} d\\\\Omega'
@Sigma_S(\\\\vec{r},\\\\hat{\\\\Omega}'\\\\cdot \\\\hat{\\\\Omega},
E_g'\\\\to E_g) @psi(\\\\vec{r},\\\\hat{\\\\Omega}',E_{g'}) +
@frac{\\\\chi_g}{4\\\\pi k} \\\\sum^G_{g'} \\\\int_{4\\\\pi}
d\\\\Omega' @nu \\\\Sigma_f(\\\\vec{r}, E_{g'}) \\\\psi(\\\\vec{r},
\\\\hat{\\\\Omega}', E_{g'}) \\\\]

where $ k $ is the eigenvalue, which physically represents the ratio
of the number of neutrons in successive generations. Operator Form In
operator form, the eigenvalue problem is

\\\\[ (\\\\mathbf{I} - \\\\mathbf{DL}^{-1}\\\\mathbf{MS})\\\\phi =
@frac{1}{k} \\\\mathbf{DL}^{-1}\\\\mathbf{M}\\\\chi \\\\mathbf{F}^T
\\\\phi \\\\, . \\\\]

where the eigenvector consists of the angular flux moments. In the
more standard form $ \\\\mathbf{A}x = k x $, we have two options. We
can keep the problem in terms of the multigroup fluxes so that

\\\\[ @overbrace{(\\\\mathbf{I} -
\\\\mathbf{DL}^{-1}\\\\mathbf{MS})^{-1}
@mathbf{DL}^{-1}\\\\mathbf{M}\\\\chi \\\\mathbf{F}^T}^{\\\\mathbf{A}}
\\\\phi = k \\\\phi \\\\, . \\\\]

Alternatively, we can solve in terms of the fission density, \\\\[ d =
\\\\mathbf{F}^T \\\\phi = \\\\sum^G_g \\\\nu \\\\Sigma_{fg} \\\\phi_g
\\\\] so that \\\\[ @overbrace{(\\\\mathbf{I} -
\\\\mathbf{DL}^{-1}\\\\mathbf{MS})^{-1}
@mathbf{DL}^{-1}\\\\mathbf{M}\\\\chi }^{\\\\mathbf{A}} d = k d \\\\, .
\\\\]

We choose the second approach, since it represents a smaller system
independent of energy and represents the more canonical formulation of
the problem. The former approach would be of value when parallelizing
over energy, as in Denovo, but that is not our focus here. A Note
About the Multigroup Solve Note that in both cases, we require the
action $ (\\\\mathbf{I} - \\\\mathbf{DL}^{-1} )$. This is equivalent
to an exact multigroup solve.

For the power method (see PowerIteration), an exact inversion is not
really necessary, since the iteration will eventually converge (though
in many more power iterations than if an exact multigroup solve were
used). In many cases, this is actually more efficient, as the number
of sweeps to get within some distance to the true eigenvector and
eigenvalue is less than if one converges on the multigroup problem
(and the within-group problems).

The only problem with this approach is that the convergence criterion
is not very well defined, since the iteration is not truly a power
iteration when exact inversions are not used (which actually means we
never really use power iteration in practice). Hence, there is no real
way to map the difference between successive iterates in the
approximate case to what the difference between successive iterates
means mathematically when full inversions are used, i.e. when a full
power iteration is used. Thus, when one does not converge the
multigroup problem, care must be taken to ensure the convergence
criterion at least somewhat represents the desired outcome.

For Krylov methods such as the Arnoldi method, stricter convergence on
the multigroup problem is required. This is because inexact (or simply
no) convergence on the multigroup problem represents an entirely
different operator $ \\\\mathbf{A} $, and hence will yield an entirely
different spectrum. Eigensolver Input Entries  \\\\[
@begin{tabular}{llll} Key & Type & Description & Default \\\\\\\\
@hline eigen\\\\_max\\\\_iters & int & Maximum iterations allowed &
100 @\\\\ eigen\\\\_tolerance & dbl & Tolerance on residual & 1e-5
@\\\\ eigen\\\\_print\\\\_out & int & Diagnostic print level & 0 @\\\\
eigen\\\\_print\\\\_interval & int & Iteration interval for
diagonostics & 10 @\\\\ @end{tabular} \\\\]

C++ includes: Eigensolver.hh ";

%feature("docstring")  detran::Eigensolver::Eigensolver "detran::Eigensolver< D >::Eigensolver(SP_mg_solver mg_solver)

Constructor.

Parameters:
-----------

mg_solver:  Multigroup solver ";

%feature("docstring")  detran::Eigensolver::~Eigensolver "detran::Eigensolver< D >::~Eigensolver()=0

Virtual destructor. ";

%feature("docstring")  detran::Eigensolver::solve "virtual void
detran::Eigensolver< D >::solve()=0

Solve the eigenvalue problem. ";


// File: classcallow_1_1EigenSolverCreator.xml
%feature("docstring") callow::EigenSolverCreator "

Creates an eigensolver.

C++ includes: EigenSolverCreator.hh ";


// File: classdetran_1_1EigenvalueManager.xml
%feature("docstring") detran::EigenvalueManager "

Manage solution of a multigroup eigenvalue problem.

C++ includes: EigenvalueManager.hh ";

%feature("docstring")  detran::EigenvalueManager::EigenvalueManager "home robertsj Research detran source src solvers EigenvalueManager cc
detran::EigenvalueManager< D >::EigenvalueManager(int argc, char
*argv[], SP_input input, SP_material material, SP_mesh mesh)

Constructor.

Parameters:
-----------

argc:  command line count

argv:  command line values

input:  parameter database

material:  material database

mesh:  mesh definition

Create the fixed source manager ";

%feature("docstring")  detran::EigenvalueManager::EigenvalueManager "detran::EigenvalueManager< D >::EigenvalueManager(SP_input input,
SP_material material, SP_mesh mesh)

Constructor (without command line)

Create the fixed source manager ";

%feature("docstring")  detran::EigenvalueManager::~EigenvalueManager "virtual detran::EigenvalueManager< D >::~EigenvalueManager()

Virtual destructor. ";

%feature("docstring")  detran::EigenvalueManager::solve "bool
detran::EigenvalueManager< D >::solve()

Solve the system. ";


// File: classdetran_1_1EnergyIndependentEigenOperator.xml
%feature("docstring") detran::EnergyIndependentEigenOperator "

Energy-independent operator for eigenvalue problems.

This operator represents the action of \\\\[ \\\\mathbf{f}^T
(\\\\mathbf{I} - \\\\mathbf{TMS})^{-1} \\\\mathbf{TM}
\\\\boldsymbol{\\\\chi} \\\\, . \\\\]

This is the preferred operator for transport eigenvalue solves since
the system is reduced in size.

C++ includes: EnergyIndependentEigenOperator.hh ";

%feature("docstring")
detran::EnergyIndependentEigenOperator::EnergyIndependentEigenOperator
"detran::EnergyIndependentEigenOperator< D
>::EnergyIndependentEigenOperator(SP_mg_solver mg_solver)

Constructor.

Parameters:
-----------

mg_solver:  Multigroup solver ";

%feature("docstring")  detran::EnergyIndependentEigenOperator::~EnergyIndependentEigenOperator 
" virtual detran::EnergyIndependentEigenOperator< D >::~EnergyIndependentEigenOperator() ";

%feature("docstring")  detran::EnergyIndependentEigenOperator::display
"void detran::EnergyIndependentEigenOperator< D >::display() const ";

%feature("docstring")
detran::EnergyIndependentEigenOperator::multiply "void
detran::EnergyIndependentEigenOperator< D >::multiply(const Vector &x,
Vector &y) ";

%feature("docstring")
detran::EnergyIndependentEigenOperator::multiply_transpose "virtual
void detran::EnergyIndependentEigenOperator< D
>::multiply_transpose(const Vector &x, Vector &y) ";


// File: classdetran_1_1Equation.xml
%feature("docstring") detran::Equation "

Traits for defining the face flux type for a discretization.

Discrete ordinates equation base.

Parameters:
-----------

D:  Problem dimension

C++ includes: Equation.hh ";

%feature("docstring")  detran::Equation::Equation "detran::Equation<
D >::Equation(SP_mesh mesh, SP_material material, SP_quadrature
quadrature, const bool update_psi)

Constructor.

Parameters:
-----------

mesh:  Geometry

material:  Material database

quadrature:  Angular mesh

update_psi:  Flag for keeping the angular flux ";

%feature("docstring")  detran::Equation::~Equation "virtual
detran::Equation< D >::~Equation() ";

%feature("docstring")  detran::Equation::solve "virtual void
detran::Equation< D >::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)=0

Solve for the cell-center and outgoing edge fluxes.

Parameters:
-----------

i:  Cell x index

j:  Cell y index

k:  Cell z index

source:  Reference to source vector

psi_in:  Incident flux for this cell

psi_out:  Outgoing flux from this cell

phi:  Reference to flux moments for this group

psi:  Reference to angular flux for this group ";

%feature("docstring")  detran::Equation::setup_group "virtual void
detran::Equation< D >::setup_group(const size_t g)=0

Setup the equations for a group.

Parameters:
-----------

g:  Current group. ";

%feature("docstring")  detran::Equation::setup_octant "virtual void
detran::Equation< D >::setup_octant(const size_t octant)=0

Setup the equations for an octant.

Parameters:
-----------

octant:  Current octant. ";

%feature("docstring")  detran::Equation::setup_angle "virtual void
detran::Equation< D >::setup_angle(const size_t angle)=0

Setup the equations for an angle.

Parameters:
-----------

angle:  Angle index within octant ";


// File: classdetran_1_1Equation__DD__1D.xml
%feature("docstring") detran::Equation_DD_1D "

Diamond difference discretization in one dimension.

See Equation_DD_3D for a general description of the diamond difference
approximation.

C++ includes: Equation_DD_1D.hh ";

%feature("docstring")  detran::Equation_DD_1D::Equation_DD_1D "home
robertsj Research detran source src transport Equation_DD_1D cc
detran::Equation_DD_1D::Equation_DD_1D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_DD_1D::solve "void
detran::Equation_DD_1D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_DD_1D::setup_group "void
detran::Equation_DD_1D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_DD_1D::setup_octant "void
detran::Equation_DD_1D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_DD_1D::setup_angle "void
detran::Equation_DD_1D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__DD__2D.xml
%feature("docstring") detran::Equation_DD_2D "

Diamond difference discretization in two dimensions.

See Equation_DD_3D for a general description of the diamond difference
approximation.

C++ includes: Equation_DD_2D.hh ";

%feature("docstring")  detran::Equation_DD_2D::Equation_DD_2D "home
robertsj Research detran source src transport Equation_DD_2D cc
detran::Equation_DD_2D::Equation_DD_2D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, const bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_DD_2D::solve "void
detran::Equation_DD_2D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_DD_2D::setup_group "void
detran::Equation_DD_2D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_DD_2D::setup_octant "void
detran::Equation_DD_2D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_DD_2D::setup_angle "void
detran::Equation_DD_2D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__DD__3D.xml
%feature("docstring") detran::Equation_DD_3D "

Diamond difference discretization in three dimensions.

Fill me.

C++ includes: Equation_DD_3D.hh ";

%feature("docstring")  detran::Equation_DD_3D::Equation_DD_3D "home
robertsj Research detran source src transport Equation_DD_3D cc
detran::Equation_DD_3D::Equation_DD_3D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_DD_3D::solve "void
detran::Equation_DD_3D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_DD_3D::setup_group "void
detran::Equation_DD_3D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_DD_3D::setup_octant "void
detran::Equation_DD_3D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_DD_3D::setup_angle "void
detran::Equation_DD_3D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__ES__1D.xml
%feature("docstring") detran::Equation_ES_1D "

Explicit slope discretization in one dimension.

The explicit slope discretization uses a weighted difference equation
of the form \\\\[ \\\\psi_{n,j} = \\\\frac{1 + \\\\alpha_{n,j}}{2}
\\\\psi_{n,j+\\\\frac{1}{2}} + \\\\frac{1 - \\\\alpha_{n,j}}{2}
\\\\psi_{n,j-\\\\frac{1}{2}} -
\\\\frac{\\\\alpha_{n,j}}{2\\\\Sigma_{tj}} T_{n, j} \\\\, , \\\\]
where \\\\[ T_{n, j} = \\\\frac{1}{2} \\\\left ( \\\\Sigma_{s,j}
\\\\hat{\\\\psi}_{n,j} + \\\\hat{Q}_{n, j} \\\\right ) \\\\, , \\\\]
contains slopes for the scattering and external source; that is, the
full right hand side is a linear function.

The weights $ \\\\alpha_{n,j} $ can be selected in several ways,
including that of the SC method, \\\\[ \\\\alpha_{n,j} = \\\\coth{
\\\\tau_{n, j} } - \\\\frac{1}{\\\\tau_{n,j}} \\\\, \\\\] and \\\\[
\\\\alpha_{n,j} = \\\\frac{\\\\tau_{n,j}}{\\\\theta + |\\\\tau_{n,j}|}
\\\\, , \\\\] where \\\\[ \\\\tau_{n, j} = \\\\frac{ \\\\Sigma_{t, j}
\\\\Delta_j }{2 \\\\mu_n} \\\\, . \\\\] In the latter, $ \\\\theta = 3
$ corresponds to a linear discontinuous, and $ \\\\theta = 1 $ to
lumped linear discontinuous. Neither of these is implemented.

The flux slope is approximated as

C++ includes: Equation_ES_1D.hh ";

%feature("docstring")  detran::Equation_ES_1D::Equation_ES_1D "detran::Equation_ES_1D::Equation_ES_1D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_ES_1D::solve "void
detran::Equation_ES_1D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_ES_1D::setup_group "void
detran::Equation_ES_1D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_ES_1D::setup_octant "void
detran::Equation_ES_1D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_ES_1D::setup_angle "void
detran::Equation_ES_1D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__MOC.xml
%feature("docstring") detran::Equation_MOC "

Method of characteristics equation base.

C++ includes: Equation_MOC.hh ";

%feature("docstring")  detran::Equation_MOC::Equation_MOC "detran::Equation_MOC::Equation_MOC(SP_mesh mesh, SP_material material,
SP_quadrature quadrature, bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_MOC::~Equation_MOC "virtual
detran::Equation_MOC::~Equation_MOC() ";

%feature("docstring")  detran::Equation_MOC::solve "virtual void
detran::Equation_MOC::solve(const size_t region, const double length,
moments_type &source, double &psi_in, double &psi_out, moments_type
&phi, angular_flux_type &psi)=0

Solve for the cell-center and outgoing edge fluxes.

Parameters:
-----------

region:  Flat source region (cardinal mesh index)

length:  Segment length

source:  Reference to sweep source vector for this group

psi_in:  Incident flux for this cell

psi_out:  Outgoing flux from this cell

phi:  Reference to flux moments for this group

psi:  Reference to angular flux for this group ";

%feature("docstring")  detran::Equation_MOC::setup_group "virtual
void detran::Equation_MOC::setup_group(const size_t g)=0

Setup the equations for a group.

Parameters:
-----------

g:  Current group. ";

%feature("docstring")  detran::Equation_MOC::setup_octant "virtual
void detran::Equation_MOC::setup_octant(const size_t o)=0

Setup the equations for an octant.

Parameters:
-----------

o:  Current octant index. ";

%feature("docstring")  detran::Equation_MOC::setup_azimuth "virtual
void detran::Equation_MOC::setup_azimuth(const size_t a)=0

Setup the equations for an azimuth.

Parameters:
-----------

a:  Azimuth within octant. ";

%feature("docstring")  detran::Equation_MOC::setup_polar "virtual
void detran::Equation_MOC::setup_polar(const size_t p)=0

Setup the equations for a polar angle.

Parameters:
-----------

p:  Polar index. ";


// File: classdetran_1_1Equation__SC__1D.xml
%feature("docstring") detran::Equation_SC_1D "

Step characteristic discretization in one dimension.

In 1D, the step characteristic approximation is the same for discrete
ordinates and MOC. The outgoing flux is defined \\\\[ \\\\psi_{out} =
A\\\\psi_{in} + B Q \\\\, , \\\\] and the average segment flux is
\\\\[ \\\\bar{\\\\psi} = \\\\frac{1}{l} \\\\Big ( B \\\\psi_{in} + C Q
\\\\Big ) \\\\, , \\\\] where \\\\[ A = e^{-\\\\Sigma_t \\\\tau} \\\\,
, \\\\] \\\\[ B = \\\\frac{1}{\\\\Sigma_t} ( 1- A ) \\\\, , \\\\] and
\\\\[ C = \\\\frac{l}{\\\\Sigma_t} \\\\Big( 1- \\\\frac{1-A}{\\\\tau}
\\\\Big ) \\\\, , \\\\] where $ l $ is the step length and $ \\\\tau =
\\\\Sigma_t l $ is optical path length.

See:   Equation_SC_MOC

C++ includes: Equation_SC_1D.hh ";

%feature("docstring")  detran::Equation_SC_1D::Equation_SC_1D "home
robertsj Research detran source src transport Equation_SC_1D cc
detran::Equation_SC_1D::Equation_SC_1D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_SC_1D::solve "void
detran::Equation_SC_1D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_SC_1D::setup_group "void
detran::Equation_SC_1D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_SC_1D::setup_octant "void
detran::Equation_SC_1D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_SC_1D::setup_angle "void
detran::Equation_SC_1D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__SC__2D.xml
%feature("docstring") detran::Equation_SC_2D "

Step characteristic discretization in two dimensions.

Reference: Lathrop, K. D. \"Spatial differencing of the Transport
Equation: Positivity vs. Accuracy\", J. Comp. Phys. 4, 475-498 (1969)

C++ includes: Equation_SC_2D.hh ";

%feature("docstring")  detran::Equation_SC_2D::Equation_SC_2D "home
robertsj Research detran source src transport Equation_SC_2D cc
detran::Equation_SC_2D::Equation_SC_2D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, const bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_SC_2D::solve "void
detran::Equation_SC_2D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_SC_2D::setup_group "void
detran::Equation_SC_2D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_SC_2D::setup_octant "void
detran::Equation_SC_2D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_SC_2D::setup_angle "void
detran::Equation_SC_2D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__SC__MOC.xml
%feature("docstring") detran::Equation_SC_MOC "

Step characteristic discretization for MOC.

In the method of characteristics, the flux is solved for along a track
assuming a flat source. For a given incident flux into a track
segment, we can define the outgoing segment flux \\\\[ \\\\psi_{out} =
A\\\\psi_{in} + B Q \\\\, , \\\\] and average segment flux \\\\[
\\\\bar{\\\\psi} = \\\\frac{1}{l} \\\\Big ( B \\\\psi_{in} + C Q
\\\\Big ) \\\\, , \\\\] where \\\\[ A = e^{-\\\\Sigma_t \\\\tau} \\\\,
, \\\\] \\\\[ B = \\\\frac{1}{\\\\Sigma_t} ( 1- A ) \\\\, , \\\\] and
\\\\[ C = \\\\frac{l}{\\\\Sigma_t} \\\\Big( 1- \\\\frac{1-A}{\\\\tau}
\\\\Big ) \\\\, , \\\\] where $ l $ is the segment length and $
\\\\tau = \\\\Sigma_t l $ is optical path length.

The step characteristic method is positive but only first-order
accurate in space.

Reference: A. Hebert, Applied Reactor Physics.

See:  Equation_DD_MOC

C++ includes: Equation_SC_MOC.hh ";

%feature("docstring")  detran::Equation_SC_MOC::Equation_SC_MOC "home
robertsj Research detran source src transport Equation_SC_MOC cc
detran::Equation_SC_MOC::Equation_SC_MOC(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, const bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_SC_MOC::solve "void
detran::Equation_SC_MOC::solve(const size_t region, const double
length, moments_type &source, double &psi_in, double &psi_out,
moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_SC_MOC::setup_group "void
detran::Equation_SC_MOC::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_SC_MOC::setup_octant "void
detran::Equation_SC_MOC::setup_octant(const size_t o)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_SC_MOC::setup_azimuth "void
detran::Equation_SC_MOC::setup_azimuth(const size_t a)

Setup the equations for an azimuth. ";

%feature("docstring")  detran::Equation_SC_MOC::setup_polar "void
detran::Equation_SC_MOC::setup_polar(const size_t p)

Setup the equations for a polar angle. ";


// File: classdetran_1_1Equation__SD__1D.xml
%feature("docstring") detran::Equation_SD_1D "

Step difference discretization in one dimension.

The step difference approximation defines \\\\[ \\\\psi_{i,n} =
\\\\left\\\\{ \\\\begin{array}{l l} \\\\psi_{i+1/2,n} & \\\\quad
\\\\text{if $\\\\mu_n > 0$} \\\\\\\\ \\\\psi_{i-1/2,n} & \\\\quad
\\\\text{if $\\\\mu_n < 0$} \\\\\\\\ \\\\end{array} \\\\right\\\\}
\\\\]

This is a first order method, but is strictly positive. Moreover, it
essentially defines a \"flat flux\" within the cell, and so it is
consistent for DGM.

C++ includes: Equation_SD_1D.hh ";

%feature("docstring")  detran::Equation_SD_1D::Equation_SD_1D "home
robertsj Research detran source src transport Equation_SD_1D cc
detran::Equation_SD_1D::Equation_SD_1D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_SD_1D::solve "void
detran::Equation_SD_1D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_SD_1D::setup_group "void
detran::Equation_SD_1D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_SD_1D::setup_octant "void
detran::Equation_SD_1D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_SD_1D::setup_angle "void
detran::Equation_SD_1D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1Equation__SD__2D.xml
%feature("docstring") detran::Equation_SD_2D "

Step difference discretization in two dimensions.

C++ includes: Equation_SD_2D.hh ";

%feature("docstring")  detran::Equation_SD_2D::Equation_SD_2D "home
robertsj Research detran source src transport Equation_SD_2D cc
detran::Equation_SD_2D::Equation_SD_2D(SP_mesh mesh, SP_material
material, SP_quadrature quadrature, const bool update_psi)

Constructor. ";

%feature("docstring")  detran::Equation_SD_2D::solve "void
detran::Equation_SD_2D::solve(const size_t i, const size_t j, const
size_t k, moments_type &source, face_flux_type &psi_in, face_flux_type
&psi_out, moments_type &phi, angular_flux_type &psi)

Solve for the cell-center and outgoing edge fluxes. ";

%feature("docstring")  detran::Equation_SD_2D::setup_group "void
detran::Equation_SD_2D::setup_group(const size_t g)

Setup the equations for a group. ";

%feature("docstring")  detran::Equation_SD_2D::setup_octant "void
detran::Equation_SD_2D::setup_octant(const size_t octant)

Setup the equations for an octant. ";

%feature("docstring")  detran::Equation_SD_2D::setup_angle "void
detran::Equation_SD_2D::setup_angle(const size_t angle)

Setup the equations for an angle. ";


// File: classdetran_1_1EquationTraits.xml
%feature("docstring") detran::EquationTraits "C++ includes:
Equation.hh ";


// File: classdetran_1_1EquationTraits_3_01__1D_01_4.xml
%feature("docstring") detran::EquationTraits< _1D > " C++ includes:
Equation.hh ";


// File: classdetran_1_1Execute.xml
%feature("docstring") detran::Execute "

Setup and execute the problem.

C++ includes: Execute.hh ";

%feature("docstring")  detran::Execute::Execute "detran::Execute::Execute(StupidParser &parser) ";

%feature("docstring")  detran::Execute::solve "template void
detran::Execute::solve< _3D >()

Solve the problem. ";

%feature("docstring")  detran::Execute::output "void
detran::Execute::output()

Write output to file. ";

%feature("docstring")  detran::Execute::dimension "int
detran::Execute::dimension() ";


// File: classdetran__external__source_1_1ExternalSource.xml
%feature("docstring") detran_external_source::ExternalSource "

Base volume source class.

Consider the general transport equation in operator form: \\\\[
\\\\mathcal{L}_g \\\\psi(\\\\vec{r}, \\\\hat{\\\\Omega}, E_g) =
Q(\\\\vec{r}, \\\\hat{\\\\Omega}, E_g) \\\\, . \\\\] The source $ Q $
is a function of space, angle, and energy. Because sources will be
used in different ways throughout Detran, a useful interface allows a
client to get discrete angular sources (for use in SN or MOC
calculations) and moment sources (for use in diffusion). For the
latter, only the isotropic component is available.

Source representation in moment form is limited to isotropic

C++ includes: ExternalSource.hh ";

%feature("docstring")
detran_external_source::ExternalSource::ExternalSource "detran_external_source::ExternalSource::ExternalSource(size_t
number_groups, SP_mesh mesh, SP_quadrature quadrature, bool
discrete=false)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

quadrature:  Pointer to angular quadrature

discrete:  Flag for discrete sources ";

%feature("docstring")
detran_external_source::ExternalSource::~ExternalSource "virtual
detran_external_source::ExternalSource::~ExternalSource()

Virtual destructor. ";

%feature("docstring")
detran_external_source::ExternalSource::number_groups "size_t
detran_external_source::ExternalSource::number_groups() const ";

%feature("docstring")
detran_external_source::ExternalSource::is_discrete "bool
detran_external_source::ExternalSource::is_discrete() const ";

%feature("docstring")  detran_external_source::ExternalSource::source
"virtual double detran_external_source::ExternalSource::source(const
size_t cell, const size_t group)=0

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran_external_source::ExternalSource::source
"virtual double detran_external_source::ExternalSource::source(const
size_t cell, const size_t group, const size_t angle)=0

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";


// File: classdetran_1_1FissionSource.xml
%feature("docstring") detran::FissionSource "

Defines the isotropic source from fission reactions.

C++ includes: FissionSource.hh ";

%feature("docstring")  detran::FissionSource::FissionSource "detran::FissionSource::FissionSource(SP_state state, SP_mesh mesh,
SP_material material)

Constructor.

Parameters:
-----------

state:   State vector

mesh:  Cartesian mesh

material:  Materials ";

%feature("docstring")  detran::FissionSource::initialize "void
detran::FissionSource::initialize()

Methods to treat fission in an outer iteration (e.g. power iteration)

Initialize to sum of cell nu*fission cross-section, normalized. ";

%feature("docstring")  detran::FissionSource::update "void
detran::FissionSource::update()

Update the fission density. ";

%feature("docstring")  detran::FissionSource::setup_outer "void
detran::FissionSource::setup_outer(const double scale=1.0)

Setup the fission source for an outer iteration.

This sets a new scaling factor $ C $ and precomputes the quantity $ v
= C \\\\times fd $.

Parameters:
-----------

scale:  Scaling factor (typically 1/keff) ";

%feature("docstring")  detran::FissionSource::source "const
State::moments_type & detran::FissionSource::source(const size_t g)

Return the fission source in a group.

The group fission source is just that component of the density
released in a particular group. Mathematically, this is just

\\\\[ q_{f,g} = @frac{\\\\chi_g}{4\\\\pi k} \\\\sum_g
\\\\nu\\\\Sigma_{f,g} \\\\phi_g \\\\, . \\\\]

Note, the scaling factor is actually arbitrary. For 2-D and 3-D, it is
$ 4\\\\pi $, possibly with the eigenvalue $ k $. The client sets this
in update.

Note also that this returns a moments source, so the client must apply
the moments-to-discrete operator.

Parameters:
-----------

g:  Group of the source.

Source vector. ";

%feature("docstring")  detran::FissionSource::density "const
State::moments_type & detran::FissionSource::density()

Return the fission density.

\\\\[ fd = \\\\sum_g \\\\nu\\\\Sigma_{f,g} \\\\phi_g \\\\, . \\\\]

Fission density vector. ";

%feature("docstring")  detran::FissionSource::set_density "void
detran::FissionSource::set_density(moments_type &f)

Set the fission density.

Parameters:
-----------

f:  User-defined density. ";

%feature("docstring")
detran::FissionSource::build_within_group_source "void
detran::FissionSource::build_within_group_source(const size_t g, const
moments_type &phi, moments_type &source)

Methods to treat fission like scatter.

Build the within group fission source.

This constructs \\\\[ q_g = \\\\chi_g \\\\nu \\\\Sigma_{fg} \\\\phi_g
\\\\, . \\\\]

Parameters:
-----------

g:  Group for this problem

phi:  Const reference to group flux.

source:  Mutable reference to moments source. ";

%feature("docstring")  detran::FissionSource::build_in_fission_source
"void detran::FissionSource::build_in_fission_source(const size_t g,
moments_type &source)

Build the in-fission source.

This constructs \\\\[ q_g = \\\\chi_g \\\\sum^G_{g',g\\\\ne g'}
\\\\Sigma_{fg'} \\\\phi_{g'} \\\\, . \\\\]

This *assumes* the state is up-to-date.

Parameters:
-----------

g:  Group for this problem

source:  Mutable reference to moments source. ";

%feature("docstring")  detran::FissionSource::build_total_group_source
"void detran::FissionSource::build_total_group_source(const size_t g,
const State::vec_moments_type &phi, State::moments_type &source)

Fill a group source vector with a fission source given a client-
defined flux vector.

For Krylov methods, we can bring the entire flux-dependent terms to
the left hand side, thus treating the flux implicitly. This function
allows the client to build the total fission source into a given group
based on an arbitrary, client-defined multigroup flux vector.

Parameters:
-----------

g:  group of source being constructed

phi:  multigroup fluxes

source:  moment vector of group source to contribute to ";

%feature("docstring")  detran::FissionSource::state "SP_state
detran::FissionSource::state()

Get the state. ";


// File: classdetran_1_1FixedSourceManager.xml
%feature("docstring") detran::FixedSourceManager "

Manage solution of a multigroup fixed source problem.

Fixed source problems are classified generically as \"fixed\" or
\"multiplying\". For the latter, the fission source is included based
on an assumed k-eigenvalue (defaulted to unity). A multiplying problem
can be solved with fission treated like scatter, or it can be added
via iteration outside the normal multigroup solver. The latter is
useful when we want to expand solutions in fission generation series.

C++ includes: FixedSourceManager.hh ";

%feature("docstring")  detran::FixedSourceManager::FixedSourceManager
"home robertsj Research detran source src solvers FixedSourceManager
cc detran::FixedSourceManager< D >::FixedSourceManager(int argc, char
*argv[], SP_input input, SP_material material, SP_mesh mesh, bool
multiply=false, bool fission=false)

Constructor.

Parameters:
-----------

argc:  command line count

argv:  command line values

input:  parameter database

material:  material database

mesh:  mesh definition

multiply:  flag for multiplying fixed source problem

fission:  ensure a fission source is built (e.g. for eigen) ";

%feature("docstring")  detran::FixedSourceManager::FixedSourceManager
"detran::FixedSourceManager< D >::FixedSourceManager(SP_input input,
SP_material material, SP_mesh mesh, bool multiply=false, bool
fission=false)

Constructor (without command line) ";

%feature("docstring")  detran::FixedSourceManager::~FixedSourceManager
"virtual detran::FixedSourceManager< D >::~FixedSourceManager()

Virtual destructor. ";

%feature("docstring")  detran::FixedSourceManager::setup "void
detran::FixedSourceManager< D >::setup()

Sets up a problem to be solved.

This sets up the discretization, quadrature, boundary, and state. By
changing the appropriate database parameters, a call to setup will
rebuild the problem using a different discretization, etc.

Source defining external sources, the quadrature is often needed. By
calling setup first, the client can extract the quadrature for
building an external source. ";

%feature("docstring")  detran::FixedSourceManager::set_source "void
detran::FixedSourceManager< D >::set_source(SP_source q)

Add a new external source.

This is called after setup, since external sources often need the
quadrature set produced by the manager. ";

%feature("docstring")  detran::FixedSourceManager::set_solver "bool
detran::FixedSourceManager< D >::set_solver()

Set the solver based on parameter database.

This must be called after setup is called and after all sources have
been set. By changing the appropriate parameters, a problem can be re-
run with a different solver. ";

%feature("docstring")  detran::FixedSourceManager::solve "bool
detran::FixedSourceManager< D >::solve(const double keff=1.0)

Solve the system.

Parameters:
-----------

keff:  Scaling factor for multiplying problems ";

%feature("docstring")  detran::FixedSourceManager::iterate "double
detran::FixedSourceManager< D >::iterate(const int generation)

Perform a fission iteration.

This function provides a way for a client to perform expansion of a
fixed source problem solution in terms of a series in fission
generations.

To use this function, the client sets up a problem in the typical way
by setting a source, options, etc. For each generation, the client
iterates. The function returns the norm of the difference between
fluxes of successive generations. Following each iteration, the state
is filled with the m-th iteration contribution to the flux. In other
words, the solution is represented as \\\\[ \\\\phi = \\\\phi_0 +
\\\\phi_1 + \\\\ldots \\\\] where $ \\\\phi_0 $ is the flux due to the
initial source, $ \\\\phi_1 $ is the \"once-fissioned\" flux, and so
on. Thus, for arbitrary eigenvalue $ k $, we have \\\\[ \\\\phi =
\\\\phi_0 + \\\\frac{1}{k}\\\\phi_1 + \\\\frac{1}{k^2}\\\\phi_2 +
\\\\ldots \\\\] Moreover, it turns out that the ratio of successive
terms tends toward $ k $ of the domain in a vacuum, $ k_v $, so that
the series can be approximately summed as \\\\[ \\\\phi = \\\\phi_0 +
\\\\frac{1}{k}\\\\phi_1 + \\\\frac{1}{k^2(1-k_v/k)}\\\\phi_2 \\\\]

Note, after generation zero, the fixed sources are eliminated. ";

%feature("docstring")  detran::FixedSourceManager::update "void
detran::FixedSourceManager< D >::update()

Update operators, etc. ";


// File: classdetran__angle_1_1GaussChebyshev.xml
%feature("docstring") detran_angle::GaussChebyshev "

Implements Gauss-Chebyshev quadrature.

The Gauss-Chebyshev quadrature approximates \\\\[ \\\\int^{1}_{-1}
\\\\frac{f(x)}{\\\\sqrt{1-x^2}} dx \\\\approx \\\\sum^m_{i=1} W_i
f(x_i) \\\\, . \\\\] Now, the integral we want is not actually
weighted, and so we set the true weights to be \\\\[ w_i = W_i
\\\\sqrt{1-x_i^2} \\\\, . \\\\]

The abscissa of the m-point quadrature are the zeros of the Chebyshev
polynomial of degree m + 1.

Relevant database parameters: quad_number_polar_octant -- number of
abscissa per half space

gausschebyshev_normalize -- normalize half space weight to 1 [false]

C++ includes: GaussChebyshev.hh ";

%feature("docstring")  detran_angle::GaussChebyshev::GaussChebyshev "home robertsj Research detran source src angle GaussChebyshev cc
detran_angle::GaussChebyshev::GaussChebyshev(const size_t
number_polar_octant, const bool normalize=false)

Constructor.

Parameters:
-----------

number_polar_octant:  number of angles per half space

normalize:  normalize half space weights to 1 ";


// File: classdetran__angle_1_1GaussLegendre.xml
%feature("docstring") detran_angle::GaussLegendre "

Implements Gauss-Legendre quadrature.

Gauss-Legendre quadrature approximates \\\\[ \\\\int^{1}_{-1} f(x) dx
\\\\approx \\\\sum^m_{i=1} w_i f(x_i) \\\\, . \\\\] The abscissa $ x_i
$ are the zeros of the Legendre polynomial of degree $ m + 1 $. The
weights are then computed numerically. Notably, an m-point G-L
quadrature exactly integrates polynomials of degree less than or equal
to 2m - 1. In general, this is the best that can be achieved when
integrating polynomial functions.

Relevant database parameters: quad_number_polar_octant -- number of
abscissa per half space

C++ includes: GaussLegendre.hh ";

%feature("docstring")  detran_angle::GaussLegendre::GaussLegendre "home robertsj Research detran source src angle GaussLegendre cc
detran_angle::GaussLegendre::GaussLegendre(const size_t
number_polar_octant)

Constructor.

Parameters:
-----------

number_polar_octant:  Number of polar angles per octant ";


// File: classcallow_1_1GaussSeidel.xml
%feature("docstring") callow::GaussSeidel "

Uses Gauss-Seidel iteration to solve a system.

Gauss-Seidel iteration in matrix form uses a splitting of the form
\\\\[ \\\\mathbf{A} = \\\\mathbf{L} + \\\\mathbf{U} + \\\\mathbf{D}
\\\\, , \\\\] which are strictly lower and upper triangle and
diagonal, respectively. The Gauss-Seidel iteration is then \\\\[
\\\\mathbf{D + L} x^{n+1} = -\\\\mathbf{U}x^{n} + b \\\\] or \\\\[
x^{n+1} = \\\\overbrace{-\\\\mathbf{D+L}^{-1}(\\\\mathbf{U})}^
{\\\\mathbf{M}}x^{n} + \\\\mathbf{D+L}^{-1}b \\\\, . \\\\]
Alternatively, one can swap $ U $ and $ L $ to produce the backward
Gauss-Seidel iteration. If used together, one has symmetric Gauss-
Seidel iteration.

It can be shown that Gauss-Seidel converges if Jacobi converges.

Because we use a sparse matrix, we actually access the elements
directly rather than via indexing.

C++ includes: GaussSeidel.hh ";

%feature("docstring")  callow::GaussSeidel::GaussSeidel "home
robertsj Research detran source src callow solver GaussSeidel cc
callow::GaussSeidel::GaussSeidel(const double atol, const double rtol,
const int maxit, const double omega=1.0, bool successive_norm=false)
";

%feature("docstring")  callow::GaussSeidel::~GaussSeidel "virtual
callow::GaussSeidel::~GaussSeidel() ";


// File: classGaussSeidelMG.xml
%feature("docstring") GaussSeidelMG "

Solves the multigroup transport equation via Gauss-Seidel.

Relevant db entries: outer_norm_type (str) [default = \"Linf\"]

C++ includes: MGSolverGS.hh ";


// File: classdetran__utilities_1_1GenException.xml
%feature("docstring") detran_utilities::GenException "

A generic mechanism to manually manage exceptions

C++ includes: GenException.hh ";

%feature("docstring")  detran_utilities::GenException::GenException "detran_utilities::GenException::GenException()

Constructs a new GenException with the default message. ";

%feature("docstring")  detran_utilities::GenException::GenException "detran_utilities::GenException::GenException(int line, std::string
file, std::string msg)

Constructs a new GenException with a provided message.

Parameters:
-----------

line:  line of code erring

file:  file in which error occurs

msg:  the message ";

%feature("docstring")  detran_utilities::GenException::what "const
char * detran_utilities::GenException::what() const  throw () Returns
the error message associated with this GenException.

the message ";

%feature("docstring")  detran_utilities::GenException::~GenException "detran_utilities::GenException::~GenException()  throw () Destroys
this GenException. ";


// File: classcallow_1_1GMRES.xml
%feature("docstring") callow::GMRES "

Uses preconditioned GMRES(m) iteration to solve a system.

GMRES seeks to find the best solution $ x $ to the linear system \\\\[
\\\\mathbf{A}x = b \\\\] such that $ x \\\\in \\\\mathcal{K}_n $,
where the Krylov subspace is defined \\\\[ \\\\mathcal{K}_n \\\\equiv
[b, \\\\mathbf{A}b, \\\\mathbf{A}^2 b, \\\\ldots, \\\\mathbf{A}^{n-1}
b] \\\\, . \\\\] Specifically, GMRES finds $ x_n $ that satisfies
\\\\[ \\\\min_{x_n \\\\in \\\\mathcal{K}_n} ||b-\\\\mathbf{A}x_n||_2
\\\\, , \\\\] i.e. it finds the least squares fit within the current
Krylov subspace at every iteration, stopping when that residual is
small enough.

In this implementation, we employ GMRES(m) as described in Kelley's
red book. A key feature is its use of Givens rotation for incremental
conversion of the upper Hessenberg matrix $ H $ to an upper triangle
matrix $ R $.

C++ includes: GMRES.hh ";

%feature("docstring")  callow::GMRES::GMRES "home robertsj Research
detran source src callow solver GMRES cc callow::GMRES::GMRES(const
double atol, const double rtol, const int maxit, const int restart=20)
";

%feature("docstring")  callow::GMRES::~GMRES "callow::GMRES::~GMRES()
";


// File: classdetran__ioutils_1_1HDF5__FileType.xml
%feature("docstring") detran_ioutils::HDF5_FileType "C++ includes:
IO_HDF5_Traits.hh ";

%feature("docstring")  detran_ioutils::HDF5_FileType::~HDF5_FileType "detran_ioutils::HDF5_FileType::~HDF5_FileType()

Destructor. Note, the type must be closed to avoid a leak. ";

%feature("docstring")  detran_ioutils::HDF5_FileType::type "hid_t
detran_ioutils::HDF5_FileType::type< double >()

Return the type for the map entry value. ";


// File: classdetran__ioutils_1_1HDF5__MemoryType.xml
%feature("docstring") detran_ioutils::HDF5_MemoryType "

Create the memory type for the map entry value as stored in memory.

Create the memory type for the map entry value as stored on disk.

It would be nice to use traits for this, but HDF5 types can be
dynamic, and so we go with a persistent container to hold that type,
which also gives a clean way to delete it.

C++ includes: IO_HDF5_Traits.hh ";

%feature("docstring")
detran_ioutils::HDF5_MemoryType::~HDF5_MemoryType "detran_ioutils::HDF5_MemoryType::~HDF5_MemoryType()

Destructor. Note, the type must be closed to avoid a leak. ";

%feature("docstring")  detran_ioutils::HDF5_MemoryType::type "hid_t
detran_ioutils::HDF5_MemoryType::type< double >()

Return the type for the map entry value. ";


// File: classMGCMTSA_1_1hh.xml
%feature("docstring") MGCMTSA::hh "

Multigroup coarse mesh transport synthetic acceleration.

C++ includes: MGCMTSA.hh ";


// File: classMGDSA_1_1hh.xml
%feature("docstring") MGDSA::hh "

Multigroup diffusion synthetic acceleration.

The multigroup DSA preconditioning process $ {P}^{-1} $ is defined to
be \\\\[ (\\\\mathbf{I} - \\\\mathbf{C}^{-1} \\\\mathbf{S}) \\\\, ,
\\\\] where $ \\\\mathbf{C} $ is the multigroup diffusion operator.
This operator treats group-to- group scattering and, if requested,
fission implicitly.

Because this performs a diffusion solve on the same mesh as the
transport problem, the resulting system can be very large. Hence, it
is likely to perform best for relatively small systems. A coarse mesh
version is under development, but this implementation will server as
the upper bound for the efficacy of diffusion-based multigroup
preconditioning.

This inherits from the shell matrix (for now) so that the action can
be used to construct an explicit operator for detailed numerical
studies

C++ includes: MGDSA.hh ";


// File: classWGDiffusionLossOperator_1_1hh.xml
%feature("docstring") WGDiffusionLossOperator::hh "

Loss operator for a one group diffusion equation.

This operator can be used for group sweeping in multigroup Krylov
diffusion solvers. Additionally, it can be used to precondition
within-group transport solves.

The one group diffusion equation is \\\\[ -\\\\nabla D(\\\\vec{r})
\\\\nabla \\\\phi(\\\\vec{r}) + \\\\Sigma_r(\\\\vec{r})
\\\\phi(\\\\vec{r}) = Q(\\\\vec{r}) \\\\] where the loss operator is
defined \\\\[ \\\\mathbf{M}[\\\\cdot] \\\\equiv (-\\\\nabla
D(\\\\vec{r}) \\\\nabla + \\\\Sigma_r(\\\\vec{r}))[\\\\cdot]\\\\, .
\\\\]

A mesh-centered discretization is used.

C++ includes: WGDiffusionLossOperator.hh ";


// File: classdetran__utilities_1_1InputDB.xml
%feature("docstring") detran_utilities::InputDB "

Flexible storage for user input.

User input for transport codes typically involves several integer and
scalar quantities along with meshing, materials data, and materials
placement. For the former quantities, it can be a bonafide pain in the
arse to maintain an input structure throughout development. To avoid
that issue, we use C++ maps of type int and double to store anything
the user needs to put there for use anywhere in the code. This is
inspired by the way Denovo handles input. For us, a descendent of
Parser will read some form of user input (be it XML or SILO or HDF5)
and fill an InputDB along with generating other things (like the data
needed to make meshes and materials libraries). Of course, some
parameters need to be there, and that can be checked during a
verification.

(It's anticipated this will become a serment-wide approach)

C++ includes: InputDB.hh ";

%feature("docstring")  detran_utilities::InputDB::InputDB "detran_utilities::InputDB::InputDB(std::string name=\"InputDB\")

Constructor with optional name. ";

%feature("docstring")  detran_utilities::InputDB::get "vec_dbl
detran_utilities::InputDB::get< vec_dbl >(const std::string &key)
const

Return value of key.

Parameters:
-----------

key:  Name of the parameter.

value:  Reference to which parameter value is assigned.

Check whether key is found. ";

%feature("docstring")  detran_utilities::InputDB::check "bool
detran_utilities::InputDB::check(const std::string &key) const

Check if the key is in a map.

Note, if used consistently, this will prevent the same key being used
for different value types. Maybe there's a better way to do this input
structure.

Parameters:
-----------

key:  Name of the parameter.

True if exists; false otherwise. ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const T value)

Put a key and value in the database.

Parameters:
-----------

key:  Name of the parameter.

value:  Reference to which parameter value is assigned.

replace:  Can we replace a current value?

Check whether key is found. ";

%feature("docstring")  detran_utilities::InputDB::get_map "const
std::map< std::string, vec_dbl > & detran_utilities::InputDB::get_map<
vec_dbl >()

Return a map. ";

%feature("docstring")  detran_utilities::InputDB::size "int
detran_utilities::InputDB::size(int type) const

Number of entries of a certain type. ";

%feature("docstring")  detran_utilities::InputDB::display "void
detran_utilities::InputDB::display() const

Display all my contents. ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const int
value) ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const double
value) ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const vec_int
value) ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const vec_dbl
value) ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const
std::string value) ";

%feature("docstring")  detran_utilities::InputDB::put "void
detran_utilities::InputDB::put(const std::string &key, const SP_input
value) ";


// File: classdetran__ioutils_1_1IO__HDF5.xml
%feature("docstring") detran_ioutils::IO_HDF5 "

Maps an HDF5 database to an InputDB and vice versa.

A problem can be specified completely using the input, a material
definition, and the mesh definition, all of which can actually live in
the InputDB (as done in StupidParser). By mapping an HDF5 file to an
input object, we can eliminate the fickle text processing, using
Python or another interface to construct the HDF5 file. Then, the
executable detran can be used, which is really handy for profiling and
other analyses.

For now, we'll use a single group, /input. All entries will go in that
group. Later, it might be useful to add ones for the mesh and material
specification.

C++ includes: IO_HDF5.hh ";

%feature("docstring")  detran_ioutils::IO_HDF5::IO_HDF5 "detran_ioutils::IO_HDF5::IO_HDF5(std::string filename)

Constructor.

Parameters:
-----------

filename:  HDF5 filename ";

%feature("docstring")  detran_ioutils::IO_HDF5::open "void
detran_ioutils::IO_HDF5::open()

Open HDF5 file for writing. This replaces old content. ";

%feature("docstring")  detran_ioutils::IO_HDF5::write "void
detran_ioutils::IO_HDF5::write(SP_input input)

Write the input database into an HDF5 file.

Parameters:
-----------

input:  Input database to be written ";

%feature("docstring")  detran_ioutils::IO_HDF5::write "void
detran_ioutils::IO_HDF5::write(SP_material mat)

Write the material database into an HDF5 file.

Parameters:
-----------

mat:  Material database to be written ";

%feature("docstring")  detran_ioutils::IO_HDF5::write "void
detran_ioutils::IO_HDF5::write(SP_mesh mesh)

Write the material database into an HDF5 file.

Parameters:
-----------

mat:  Material database to be written ";

%feature("docstring")  detran_ioutils::IO_HDF5::close "void
detran_ioutils::IO_HDF5::close()

Close an HDF5 file if open. ";

%feature("docstring")  detran_ioutils::IO_HDF5::read_input "SP_input
detran_ioutils::IO_HDF5::read_input() ";

%feature("docstring")  detran_ioutils::IO_HDF5::read_material "SP_material detran_ioutils::IO_HDF5::read_material()

Get a material database from file. ";

%feature("docstring")  detran_ioutils::IO_HDF5::read_mesh "SP_mesh
detran_ioutils::IO_HDF5::read_mesh()

Get a mesh from file. ";


// File: classdetran__external__source_1_1IsotropicSource.xml
%feature("docstring") detran_external_source::IsotropicSource "

Isotropic volume source.

C++ includes: IsotropicSource.hh ";

%feature("docstring")
detran_external_source::IsotropicSource::IsotropicSource "home
robertsj Research detran source src external_source IsotropicSource cc
detran_external_source::IsotropicSource::IsotropicSource(size_t
number_groups, SP_mesh mesh, spectra_type &spectra, vec_int &map,
SP_quadrature quadrature=SP_quadrature(0))

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

spectra:  Vector of spectra, size = [#spectra][#groups]. The unit is n
/cc-sec.

map:  Map of where spectra are located, size = [#cells]

quadrature:  Pointer to quadrature (optional) ";

%feature("docstring")  detran_external_source::IsotropicSource::source
"double detran_external_source::IsotropicSource::source(const size_t
cell, const size_t group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran_external_source::IsotropicSource::source
"double detran_external_source::IsotropicSource::source(const size_t
cell, const size_t group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";


// File: classcallow_1_1Jacobi.xml
%feature("docstring") callow::Jacobi "

Uses Jacobi iteration to solve a system.

Jacobi iteration in matrix form uses a splitting of the form \\\\[
\\\\mathbf{A} = \\\\mathbf{L} + \\\\mathbf{U} + \\\\mathbf{D} \\\\, ,
\\\\] which are strictly lower and upper triangle and diagonal,
respectively. The Jacobi iteration is then \\\\[ \\\\mathbf{D} x^{n+1}
= -(\\\\mathbf{L}+\\\\mathbf{U})x^{n} + b \\\\] or \\\\[ x^{n+1} =
\\\\overbrace{-\\\\mathbf{D}^{-1}(\\\\mathbf{L}+\\\\mathbf{U})}^
{\\\\mathbf{M}}x^{n} + \\\\mathbf{D}^{-1}b \\\\, . \\\\] The procedure
converges if the iteration matrix is bounded, i.e. \\\\[
\\\\rho(\\\\mathbf{M}) < 1 \\\\, , \\\\] which is guaranteed if $
\\\\mathbf{A} $ is strictly diagonally dominant, meaning that \\\\[
|a_{ii}| > \\\\sum_{j \\\\neq i} |a_{ij}| \\\\, . \\\\]

Because we use a sparse matrix, we actually access the elements
directly rather than via indexing.

C++ includes: Jacobi.hh ";

%feature("docstring")  callow::Jacobi::Jacobi "home robertsj Research
detran source src callow solver Jacobi cc callow::Jacobi::Jacobi(const
double atol, const double rtol, const int maxit, bool
successive_norm=false) ";

%feature("docstring")  callow::Jacobi::~Jacobi "virtual
callow::Jacobi::~Jacobi() ";


// File: classJFNK.xml
%feature("docstring") JFNK "

Solves the eigenvalue problem via Jacobian-Free Newton-Krylov.

The eigenvalue problem can be cast in the nonlinear form

f = | (A-k*I)*fd | = 0 | 1/2 - fd'*fd |

where A is the monoenergetic eigenoperator.

The Jacobian is

J = | (A-k*I) fd | = 0 | fd' 0 |

C++ includes: JFNK.hh ";


// File: classdetran_1_1KineticsMaterial.xml
%feature("docstring") detran::KineticsMaterial "

Extends Material class for use in kinetics problems.

The basic structure of the materials class for time-dependent problems
is to expose the typical Material interface, extend it to include
basic kinetics parameters, and to provide an interface for defining
materials at a time. For all cases, base cross section data is
defined. For a few parameters, a synthetic value is defined that
accounts for the time discretization. While in the most general case,
materials will be dependent on the fine mesh, a region base material
is likely present. For the case of a fine mesh-dependent material, the
user will have to explicitly create a material map with unique id's
for each fine mesh.

Time-dependent transport in reactor applications requires tracking of
delayed neutrons. This class contains the delayed neutron fractions
(i.e. $ \\\\beta_i $) for each delayed precursor group $ i $ for every
material.

Some assumptions: assume that $ v_g $ is time-independent

assume that $ \\\\beta_i $ is material-dependent

assume that $ \\\\lambda_i $ is material-independent

C++ includes: KineticsMaterial.hh ";

%feature("docstring")  detran::KineticsMaterial::KineticsMaterial "home robertsj Research detran source src kinetics KineticsMaterial cc
detran::KineticsMaterial::KineticsMaterial(const size_t
number_materials, const size_t number_energy_groups, const size_t
number_precursor_groups, std::string name=\"KineticsMaterial\")

Constructor.

Parameters:
-----------

number_materials:  Number of materials

number_energy_groups:  Number of energy groups

number_precursor_groups:  Number of precursor groups

name:  Material name ";

%feature("docstring")  detran::KineticsMaterial::~KineticsMaterial "virtual detran::KineticsMaterial::~KineticsMaterial()

Virtual destructor. ";

%feature("docstring")  detran::KineticsMaterial::set_velocity "void
detran::KineticsMaterial::set_velocity(const size_t g, const double v)

Set the velocity for group $ i $. ";

%feature("docstring")  detran::KineticsMaterial::set_lambda "void
detran::KineticsMaterial::set_lambda(const size_t i, const double v)

Set the decay constant for group $ i $. ";

%feature("docstring")  detran::KineticsMaterial::set_beta "void
detran::KineticsMaterial::set_beta(const size_t m, const size_t i,
const double v)

Set delayed neutron fraction for a material and group. ";

%feature("docstring")  detran::KineticsMaterial::set_beta "void
detran::KineticsMaterial::set_beta(const size_t i, const double v)

Set delayed neutron fraction for a group (constant for all materials)
";

%feature("docstring")  detran::KineticsMaterial::set_chi_d "void
detran::KineticsMaterial::set_chi_d(const size_t m, const size_t i,
const size_t g, double v)

Set delayed neutron spectrum. ";

%feature("docstring")  detran::KineticsMaterial::velocity "double
detran::KineticsMaterial::velocity(const size_t g) const

Get the velocity for group $ g $. ";

%feature("docstring")  detran::KineticsMaterial::lambda "double
detran::KineticsMaterial::lambda(const size_t i) const

Get the decay constant for group $ i $. ";

%feature("docstring")  detran::KineticsMaterial::beta "double
detran::KineticsMaterial::beta(const size_t m, const size_t i) const

Get delayed neutron fraction for a material and group. ";

%feature("docstring")  detran::KineticsMaterial::beta_total "double
detran::KineticsMaterial::beta_total(const size_t m) const

Get total delayed neutron fraction for a material. ";

%feature("docstring")  detran::KineticsMaterial::chi_d "double
detran::KineticsMaterial::chi_d(const size_t m, const size_t i, const
size_t g) const

Get delayed neutron spectrum. ";

%feature("docstring")
detran::KineticsMaterial::number_precursor_groups "KineticsMaterial::size_t
detran::KineticsMaterial::number_precursor_groups() const

Get the number of precursor groups. ";

%feature("docstring")  detran::KineticsMaterial::finalize "void
detran::KineticsMaterial::finalize()

Check that all values are positive, etc. ";

%feature("docstring")  detran::KineticsMaterial::display "void
detran::KineticsMaterial::display()

Display. ";


// File: classdetran__angle_1_1LevelSymmetric.xml
%feature("docstring") detran_angle::LevelSymmetric "

2D/3D Level-symmetric (LQn) quadrature class.

Level symmetric quadratures are characterized by using the same set of
$ N/2 $ positive values of the direction cosines with respect to each
of the axes. There are $ N(N+2)/8 $ ordinates per octant, yielding $
N(N+2) $ directions for a 3D problem.

Not all the direction cosines are independent; for a given $ \\\\mu_n
$, one has a single degree of freedom, e.g. one can further choose
only $ \\\\eta_n $. The final cosine is of course defined such that
the sum in quadrature is unity.

An unfortunate aspect of LQn is that negative weights appear for $ N =
20 $. As an alternative, see the uniform and equal weight quadrature
set (UEn) UniformEqual, which yields positive weights for arbitrarily
high $ N $.

Relevant database parameters: quad_number_polar_octant -- number of
abscissa per half space

C++ includes: LevelSymmetric.hh ";

%feature("docstring")  detran_angle::LevelSymmetric::LevelSymmetric "home robertsj Research detran source src angle LevelSymmetric cc
detran_angle::LevelSymmetric::LevelSymmetric(size_t np, size_t dim)

Constructor.

Parameters:
-----------

np:  Number of abscissa per axis in one octant (2*np = order)

dim:  Problem dimension ";


// File: classdetran_1_1LinearExternalSource.xml
%feature("docstring") detran::LinearExternalSource "

Base class for time-dependent external sources.

The only addition beyond the base external source class is to specify
a time at which the source will be evaluated during the transport
solve for the current step.

C++ includes: LinearExternalSource.hh ";

%feature("docstring")
detran::LinearExternalSource::LinearExternalSource "home robertsj
Research detran source src kinetics LinearExternalSource cc
detran::LinearExternalSource::LinearExternalSource(const size_t
number_groups, SP_mesh mesh, vec_dbl times, vec_source sources, bool
discrete=false)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

times:  Times at which sources are defined

sources:  Sources at each time

discrete:  Flag to indicate treatment as moment or discrete ";

%feature("docstring")
detran::LinearExternalSource::~LinearExternalSource "virtual
detran::LinearExternalSource::~LinearExternalSource()

Virtual destructor. ";

%feature("docstring")  detran::LinearExternalSource::source "double
detran::LinearExternalSource::source(const size_t cell, const size_t
group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran::LinearExternalSource::source "double
detran::LinearExternalSource::source(const size_t cell, const size_t
group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";

%feature("docstring")  detran::LinearExternalSource::set_time "void
detran::LinearExternalSource::set_time(const double time)

Set the time of the source.

Parameters:
-----------

time:  Time at which source is evaluated ";


// File: classdetran_1_1LinearMaterial.xml
%feature("docstring") detran::LinearMaterial "

Material defined via discrete materials to be linearly interpolated.

The user specified material databases for a sequence of monotonically
increasing times, \\\\[ t_0, t_1, \\\\cdots, t_N \\\\quad t_0 \\\\ge 0
\\\\] If the initial time is not zero, then for times times $ t \\\\le
t_0 $, the first material is returned. Likewise, if $t  t_N$, the last
material is returned.

C++ includes: LinearMaterial.hh ";

%feature("docstring")  detran::LinearMaterial::LinearMaterial "home
robertsj Research detran source src kinetics LinearMaterial cc
detran::LinearMaterial::LinearMaterial(const vec_dbl &times, const
vec_material &materials, std::string name=\"LinearMaterial\")

Constructor.

Parameters:
-----------

times:  Times at which materials are defined

materials:  Materials defined at discrete times

The first time *must* be zero ";

%feature("docstring")  detran::LinearMaterial::~LinearMaterial "virtual detran::LinearMaterial::~LinearMaterial()

Virtual destructor. ";


// File: classcallow_1_1LinearSolver.xml
%feature("docstring") callow::LinearSolver "

Base class for iterative linear solvers.

We solve problems of the form \\\\[ \\\\mathbf{A}x = b \\\\] via
iterative methods. A system is \"solved\" when the norm of the
residual is small enough or some maximum iteration count is reached.
The residual is defined \\\\[ \\\\mathbf{A}x - b \\\\] and its norm is
\\\\[ r = || \\\\mathbf{A}x - b || \\\\] By default, the L2 norm is
used, though L1 and Linf are also recorded can can be used. This
represents the absolute norm. Sometimes it makes more sense to check
the residual with respect to the initial norm, in which case the
relative norm is \\\\[ r_n / r_0 = || \\\\mathbf{A}x^{n} - b || / ||
\\\\mathbf{A}x^{0} - b || \\\\] The iteration terminates when \\\\[
r_n < \\\\mathrm{max} ( \\\\tau_{\\\\mathrm{rel}} r_0,
\\\\tau_{\\\\mathrm{abs}} ) \\\\]

The solvers currently implemented are  Richardson

Jacobi

Gauss-Seidel

GMRES(m) along with Jacobi and ILU0 preconditioners. If PETSc is
enabled, all of its solvers are potentially available.

Note, some linear solvers require that the matrix provides L, U, and D
operations. Here, we simply require those solvers to have a Matrix
operator or subclasses so that the elements can be accessed directly.

C++ includes: LinearSolver.hh ";

%feature("docstring")  callow::LinearSolver::LinearSolver "callow::LinearSolver::LinearSolver(const double atol, const double
rtol, const int maxit=100, std::string name=\"solver\") ";

%feature("docstring")  callow::LinearSolver::~LinearSolver "virtual
callow::LinearSolver::~LinearSolver() ";

%feature("docstring")  callow::LinearSolver::set_operators "void
callow::LinearSolver::set_operators(SP_matrix A, SP_db db=SP_db(0))

Sets the operators for the linear system to solve.

Parameters:
-----------

A:  linear operator

P:  optional preconditioning process

side:  specifies on what side of A the preconditioner operates ";

%feature("docstring")  callow::LinearSolver::set_preconditioner "virtual void
callow::LinearSolver::set_preconditioner(SP_preconditioner P, const
int side=LEFT)

Set the preconditioner. This allows the client to build, change, etc.
";

%feature("docstring")  callow::LinearSolver::preconditioner "SP_preconditioner callow::LinearSolver::preconditioner()

Get the preconditioner. This allows the client to build, change the
curent PC.s ";

%feature("docstring")  callow::LinearSolver::set_tolerances "void
callow::LinearSolver::set_tolerances(const double atol, const double
rtol, const int maxit)

Parameters:
-----------

atol:  absolute tolerance (||r_n|| < atol)

rtol:  relative tolerance (||r_n|| < rtol * ||r_0||)

maxit:  maximum iterations (n < maxit) ";

%feature("docstring")  callow::LinearSolver::set_monitor_level "void
callow::LinearSolver::set_monitor_level(const int v)

Print residual norms and other diagonostic information.

Parameters:
-----------

v:  monitor via stdout ";

%feature("docstring")  callow::LinearSolver::set_monitor_diverge "void callow::LinearSolver::set_monitor_diverge(const bool v)

Turn on monitoring of diverging iterations. ";

%feature("docstring")  callow::LinearSolver::set_norm_type "void
callow::LinearSolver::set_norm_type(const int norm_type)

Set a norm type. ";

%feature("docstring")  callow::LinearSolver::solve "int
callow::LinearSolver::solve(const Vector &b, Vector &x)

Parameters:
-----------

b:  right hand side

x:  unknown vector ";

%feature("docstring")  callow::LinearSolver::residual_norms "std::vector<double> callow::LinearSolver::residual_norms()

return the residual norms ";

%feature("docstring")  callow::LinearSolver::number_iterations "int
callow::LinearSolver::number_iterations() const

return the number of iterations ";


// File: classcallow_1_1LinearSolverCreator.xml
%feature("docstring") callow::LinearSolverCreator "

Creates a.

C++ includes: LinearSolverCreator.hh ";


// File: classdetran__user_1_1LRA.xml
%feature("docstring") detran_user::LRA "C++ includes: LRA.hh ";

%feature("docstring")  detran_user::LRA::LRA "home robertsj Research
detran source src solvers time LRA cc detran_user::LRA::LRA(SP_mesh
mesh, bool doingtransport, bool steady)

Constructor.

Parameters:
-----------

mesh:  User-defined LRA mesh ";

%feature("docstring")  detran_user::LRA::set_state "void
detran_user::LRA::set_state(SP_state)

Set the state vector. ";

%feature("docstring")  detran_user::LRA::initialize_materials "void
detran_user::LRA::initialize_materials() ";

%feature("docstring")  detran_user::LRA::update_P_and_T "void
detran_user::LRA::update_P_and_T(double t, double dt) ";

%feature("docstring")  detran_user::LRA::T "vec_dbl
detran_user::LRA::T() ";

%feature("docstring")  detran_user::LRA::P "vec_dbl
detran_user::LRA::P() ";

%feature("docstring")  detran_user::LRA::physics "SP_multiphysics
detran_user::LRA::physics() ";

%feature("docstring")  detran_user::LRA::area "double
detran_user::LRA::area() ";

%feature("docstring")  detran_user::LRA::set_area "void
detran_user::LRA::set_area(double a) ";

%feature("docstring")  detran_user::LRA::update_impl "void
detran_user::LRA::update_impl() ";


// File: classdetran_1_1Manager.xml
%feature("docstring") detran::Manager "

Manager for initializing and finalizing external libraries.

This is useful if the user wants finer control on using detran classes
than PyExecute affords but still requires libraries to be initialized
and finalized.

C++ includes: Manager.hh ";


// File: classdetran__material_1_1Material.xml
%feature("docstring") detran_material::Material "

Simple cross section container.

All data is stored with the material index changing fastest. This
appears to be the best storage scheme with respect to memory access.

C++ includes: Material.hh ";

%feature("docstring")  detran_material::Material::Material "detran_material::Material::Material(const size_t number_materials,
const size_t number_groups, std::string name=\"no name given\")

Constructor.

Parameters:
-----------

number_materials:  Number of materials.

number_groups:  Number of energy groups.

downscatter:  Switch on to use only downscatter. ";

%feature("docstring")  detran_material::Material::~Material "virtual
detran_material::Material::~Material()

Virtual destructor. ";

%feature("docstring")  detran_material::Material::set_downscatter "void detran_material::Material::set_downscatter(bool v)

Explicitly turn on downscatter-only. ";

%feature("docstring")  detran_material::Material::set_sigma_t "void
detran_material::Material::set_sigma_t(size_t m, size_t g, double v)
";

%feature("docstring")  detran_material::Material::set_sigma_a "void
detran_material::Material::set_sigma_a(size_t m, size_t g, double v)
";

%feature("docstring")  detran_material::Material::set_nu_sigma_f "void detran_material::Material::set_nu_sigma_f(size_t m, size_t g,
double v) ";

%feature("docstring")  detran_material::Material::set_sigma_f "void
detran_material::Material::set_sigma_f(size_t m, size_t g, double v)
";

%feature("docstring")  detran_material::Material::set_nu "void
detran_material::Material::set_nu(size_t m, size_t g, double v) ";

%feature("docstring")  detran_material::Material::set_chi "void
detran_material::Material::set_chi(size_t m, size_t g, double v) ";

%feature("docstring")  detran_material::Material::set_sigma_s "void
detran_material::Material::set_sigma_s(size_t m, size_t g, size_t gp,
double v) ";

%feature("docstring")  detran_material::Material::set_diff_coef "void
detran_material::Material::set_diff_coef(size_t m, size_t g, double v)
";

%feature("docstring")  detran_material::Material::set_sigma_t "void
detran_material::Material::set_sigma_t(size_t m, vec_dbl &v) ";

%feature("docstring")  detran_material::Material::set_sigma_a "void
detran_material::Material::set_sigma_a(size_t m, vec_dbl &v) ";

%feature("docstring")  detran_material::Material::set_nu_sigma_f "void detran_material::Material::set_nu_sigma_f(size_t m, vec_dbl &v)
";

%feature("docstring")  detran_material::Material::set_sigma_f "void
detran_material::Material::set_sigma_f(size_t m, vec_dbl &v) ";

%feature("docstring")  detran_material::Material::set_nu "void
detran_material::Material::set_nu(size_t m, vec_dbl &v) ";

%feature("docstring")  detran_material::Material::set_chi "void
detran_material::Material::set_chi(size_t m, vec_dbl &v) ";

%feature("docstring")  detran_material::Material::set_sigma_s "void
detran_material::Material::set_sigma_s(size_t m, size_t g, vec_dbl &v)
";

%feature("docstring")  detran_material::Material::set_diff_coef "void
detran_material::Material::set_diff_coef(size_t m, vec_dbl &v) ";

%feature("docstring")  detran_material::Material::sigma_t "double
detran_material::Material::sigma_t(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::sigma_a "double
detran_material::Material::sigma_a(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::nu_sigma_f "double
detran_material::Material::nu_sigma_f(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::sigma_f "double
detran_material::Material::sigma_f(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::nu "double
detran_material::Material::nu(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::chi "double
detran_material::Material::chi(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::sigma_s "double
detran_material::Material::sigma_s(size_t m, size_t g, size_t gp)
const ";

%feature("docstring")  detran_material::Material::diff_coef "double
detran_material::Material::diff_coef(size_t m, size_t g) const ";

%feature("docstring")  detran_material::Material::sigma_t "Material::vec_dbl detran_material::Material::sigma_t(size_t m) const
";

%feature("docstring")  detran_material::Material::sigma_a "Material::vec_dbl detran_material::Material::sigma_a(size_t m) const
";

%feature("docstring")  detran_material::Material::nu_sigma_f "Material::vec_dbl detran_material::Material::nu_sigma_f(size_t m)
const ";

%feature("docstring")  detran_material::Material::sigma_f "Material::vec_dbl detran_material::Material::sigma_f(size_t m) const
";

%feature("docstring")  detran_material::Material::nu "Material::vec_dbl detran_material::Material::nu(size_t m) const ";

%feature("docstring")  detran_material::Material::chi "Material::vec_dbl detran_material::Material::chi(size_t m) const ";

%feature("docstring")  detran_material::Material::sigma_s "Material::vec2_dbl detran_material::Material::sigma_s(size_t m) const
";

%feature("docstring")  detran_material::Material::diff_coef "Material::vec_dbl detran_material::Material::diff_coef(size_t m) const
";

%feature("docstring")  detran_material::Material::number_groups "size_t detran_material::Material::number_groups() const ";

%feature("docstring")  detran_material::Material::number_materials "size_t detran_material::Material::number_materials() const ";

%feature("docstring")  detran_material::Material::lower "Material::size_t detran_material::Material::lower(size_t g) const

Lower scatter group bound.

This is the *lowest* index (highest energy) $ g' $ that leads to
downscatter for a given outgoing group $ g $. ";

%feature("docstring")  detran_material::Material::upper "Material::size_t detran_material::Material::upper(size_t g) const

Upper scatter group bound.

This is the *highest* index (lowest energy) $ g' $ that upscatters
size_to the outgoing group $ g $. ";

%feature("docstring")  detran_material::Material::downscatter "bool
detran_material::Material::downscatter()

Do we do only downscatter? ";

%feature("docstring")  detran_material::Material::upscatter_cutoff "size_t detran_material::Material::upscatter_cutoff()

Index below which upscatter doesn't occur for any material. ";

%feature("docstring")  detran_material::Material::compute_sigma_a "void detran_material::Material::compute_sigma_a()

Compute the absorption cross section from total and scattering.

this overwrites any data for $ \\\\Sigma_a $ already stored. ";

%feature("docstring")  detran_material::Material::compute_diff_coef "void detran_material::Material::compute_diff_coef()

Compute the diffusion coefficient from $ \\\\Sigma_t $.

Assuming isotropic scattering in the LAB, the diffusion coefficient is
simply $ D = 1/3\\\\Sigma_t $.

Todo Update diffusion definition if anisotropic scattering is added.

This overwrites any data for $ D $ already stored. ";

%feature("docstring")  detran_material::Material::finalize "void
detran_material::Material::finalize()

Computes scattering bounds and absorption cross section. ";

%feature("docstring")  detran_material::Material::display "void
detran_material::Material::display()

Pretty print the material database. ";


// File: classcallow_1_1Matrix.xml
%feature("docstring") callow::Matrix "

CRS matrix.

This is the base matrix class used within callow. It implements a
compressed row storage (CRS) matrix. An example of what this means is
as follows.

Example:

| 7 0 0 2 | | 0 2 0 4 | | 1 0 0 0 | | 3 8 0 6 |

value = [7 2 2 4 1 3 8 6] column indices = [0 3 1 3 0 0 1 3] row
pointers = [0 2 4 5 8]

We're keeping this simple to use. The use must specify the number of
nonzero entries per row, either via one value applied to all rows or
an array of values for each row. The reason the row storage is needed
rather than say the total number of nonzero entries is to make
construction easier. During construction, the matrix is stored
(temporarily) in coordinate (COO) format, i.e. a list of (i, j, value)
triplets. If this were stored in one monolithic array, the
construction of the CSR structure would require we sort all the
triplets by row and then by column (since we want explicit access to
L, D, and U). Initial testing proved that sorting is just too time
consuming, even with an n*log(n) method like quicksort.

Consequently, we do keep (i, j, value) triplets, but they are stored
by row, and to store by row, we need an initial guess of how many
entries there are. For now, the size can not be increased, but that
would not be too difficult.

During the construction process, a COO (row, column, value) format is
used. This allows the user to add entries one at a time, a row at a
time, a column at a time, or several rcv triples at a time. At the
construction process, the storage is streamlined into CSR format,
ensuring that entries are stored by row and column and that a diagonal
entry exists. That latter is required for things like the Jacobi or
GaussSeidel solvers, along with certain preconditioner types.

C++ includes: Matrix.hh ";

%feature("docstring")  callow::Matrix::Matrix "home robertsj Research
detran source src callow matrix Matrix cc callow::Matrix::Matrix() ";

%feature("docstring")  callow::Matrix::Matrix "callow::Matrix::Matrix(const int m, const int n) ";

%feature("docstring")  callow::Matrix::Matrix "callow::Matrix::Matrix(const int m, const int n, const int nnz) ";

%feature("docstring")  callow::Matrix::Matrix "callow::Matrix::Matrix(Matrix &A) ";

%feature("docstring")  callow::Matrix::~Matrix "callow::Matrix::~Matrix() ";

%feature("docstring")  callow::Matrix::preallocate "void
callow::Matrix::preallocate(const int nnz_row)

allocate using constant row size ";

%feature("docstring")  callow::Matrix::preallocate "void
callow::Matrix::preallocate(int *nnz_rows)

allocate using variable row size ";

%feature("docstring")  callow::Matrix::insert "bool
callow::Matrix::insert(int i, int j, double v, const int type=INSERT)

add one value (return false if can't add) ";

%feature("docstring")  callow::Matrix::insert "bool
callow::Matrix::insert(int i, int *j, double *v, int n, const int
type=INSERT)

add n values to a row (return false if can't add) ";

%feature("docstring")  callow::Matrix::insert "bool
callow::Matrix::insert(int *i, int j, double *v, int n, const int
type=INSERT)

add n values to a column (return false if can't add) ";

%feature("docstring")  callow::Matrix::insert "bool
callow::Matrix::insert(int *i, int *j, double *v, int n, const int
type=INSERT)

add n triplets (return false if can't add) ";

%feature("docstring")  callow::Matrix::start "int
callow::Matrix::start(const int i) const

starting index for a row ";

%feature("docstring")  callow::Matrix::diagonal "int
callow::Matrix::diagonal(const int i) const

diagonal index for a row ";

%feature("docstring")  callow::Matrix::end "int
callow::Matrix::end(const int i) const

ending index for a row ";

%feature("docstring")  callow::Matrix::column "int
callow::Matrix::column(const int p) const

column index from cardinal index ";

%feature("docstring")  callow::Matrix::values "double*
callow::Matrix::values() ";

%feature("docstring")  callow::Matrix::columns "int*
callow::Matrix::columns() ";

%feature("docstring")  callow::Matrix::rows "int*
callow::Matrix::rows() ";

%feature("docstring")  callow::Matrix::diagonals "int*
callow::Matrix::diagonals() ";

%feature("docstring")  callow::Matrix::number_nonzeros "int
callow::Matrix::number_nonzeros() const

number of nonzeros ";

%feature("docstring")  callow::Matrix::allocated "bool
callow::Matrix::allocated() const

is memory allocated? ";

%feature("docstring")  callow::Matrix::print_matlab "void
callow::Matrix::print_matlab(std::string filename=\"matrix.out\")
const

print (i, j, v) to ascii file with 1-based indexing for matlab ";

%feature("docstring")  callow::Matrix::assemble "void
callow::Matrix::assemble() ";

%feature("docstring")  callow::Matrix::multiply "void
callow::Matrix::multiply(const Vector &x, Vector &y) ";

%feature("docstring")  callow::Matrix::multiply_transpose "void
callow::Matrix::multiply_transpose(const Vector &x, Vector &y) ";

%feature("docstring")  callow::Matrix::display "void
callow::Matrix::display() const ";


// File: classcallow_1_1MatrixBase.xml
%feature("docstring") callow::MatrixBase "C++ includes: MatrixBase.hh
";

%feature("docstring")  callow::MatrixBase::MatrixBase "callow::MatrixBase::MatrixBase() ";

%feature("docstring")  callow::MatrixBase::MatrixBase "callow::MatrixBase::MatrixBase(const int m, const int n) ";

%feature("docstring")  callow::MatrixBase::~MatrixBase "virtual
callow::MatrixBase::~MatrixBase() ";

%feature("docstring")  callow::MatrixBase::set_size "virtual void
callow::MatrixBase::set_size(const int m, const int n) ";

%feature("docstring")  callow::MatrixBase::number_rows "int
callow::MatrixBase::number_rows() const ";

%feature("docstring")  callow::MatrixBase::number_columns "int
callow::MatrixBase::number_columns() const ";

%feature("docstring")  callow::MatrixBase::is_ready "bool
callow::MatrixBase::is_ready() const ";

%feature("docstring")  callow::MatrixBase::multiply "void
callow::MatrixBase::multiply(SP_vector x, SP_vector y) ";

%feature("docstring")  callow::MatrixBase::multiply_transpose "void
callow::MatrixBase::multiply_transpose(SP_vector x, SP_vector y) ";

%feature("docstring")  callow::MatrixBase::compute_explicit "virtual
void callow::MatrixBase::compute_explicit(std::string
filename=\"matrix.out\") ";

%feature("docstring")  callow::MatrixBase::assemble "virtual void
callow::MatrixBase::assemble()=0 ";

%feature("docstring")  callow::MatrixBase::multiply "virtual void
callow::MatrixBase::multiply(const Vector &x, Vector &y)=0 ";

%feature("docstring")  callow::MatrixBase::multiply_transpose "virtual void callow::MatrixBase::multiply_transpose(const Vector &x,
Vector &y)=0 ";

%feature("docstring")  callow::MatrixBase::display "virtual void
callow::MatrixBase::display() const =0 ";

%feature("docstring")  callow::MatrixBase::print_matlab "virtual void
callow::MatrixBase::print_matlab(std::string filename=\"matrix.out\")
const ";


// File: classcallow_1_1MatrixDense.xml
%feature("docstring") callow::MatrixDense "

Dense matrix.

This is a simple matrix stored in row-major format. It should be easy
enough to use BLAS routines.

Note, PETSc is accessible, but note that PETSc uses a column-major
order (like Fortran). Hence, we switch the sizes and use transpose by
default.

C++ includes: MatrixDense.hh ";

%feature("docstring")  callow::MatrixDense::MatrixDense "home
robertsj Research detran source src callow matrix MatrixDense cc
callow::MatrixDense::MatrixDense(const int m, const int n, const
double v=0.0) ";

%feature("docstring")  callow::MatrixDense::MatrixDense "callow::MatrixDense::MatrixDense(const MatrixDense &A) ";

%feature("docstring")  callow::MatrixDense::~MatrixDense "callow::MatrixDense::~MatrixDense() ";

%feature("docstring")  callow::MatrixDense::insert "bool
callow::MatrixDense::insert(int i, int j, double v, const int
type=INSERT)

add one value (return false if can't add) ";

%feature("docstring")  callow::MatrixDense::insert_row "bool
callow::MatrixDense::insert_row(int i, double *v, const int
type=INSERT)

add a row (return false if can't add) ";

%feature("docstring")  callow::MatrixDense::insert_col "bool
callow::MatrixDense::insert_col(int j, double *v, const int
type=INSERT)

add a column (return false if can't add) ";

%feature("docstring")  callow::MatrixDense::values "double*
callow::MatrixDense::values() ";

%feature("docstring")  callow::MatrixDense::print_matlab "void
callow::MatrixDense::print_matlab(std::string filename=\"matrix.out\")
const

print (i, j, v) to ascii file with 1-based indexing for matlab ";

%feature("docstring")  callow::MatrixDense::assemble "void
callow::MatrixDense::assemble() ";

%feature("docstring")  callow::MatrixDense::multiply "void
callow::MatrixDense::multiply(const Vector &x, Vector &y) ";

%feature("docstring")  callow::MatrixDense::multiply_transpose "void
callow::MatrixDense::multiply_transpose(const Vector &x, Vector &y) ";

%feature("docstring")  callow::MatrixDense::display "void
callow::MatrixDense::display() const ";


// File: classcallow_1_1MatrixShell.xml
%feature("docstring") callow::MatrixShell "

Defines a matrix free operator.

For many iterative methods, only the action of the operator on a
vector is required. Frequently, constructing a matrix explicitly is
too memory intensive, and so a matrix free operator that defines the
action is desirable.

This class is abstract. The client must define the action of the
matrix and its transpose (though for the latter, the client may simply
throw an exception to forbid its use).

C++ includes: MatrixShell.hh ";

%feature("docstring")  callow::MatrixShell::MatrixShell "home
robertsj Research detran source src callow matrix MatrixShell cc
callow::MatrixShell::MatrixShell(void *context)

the context is the \"this\" of the caller ";

%feature("docstring")  callow::MatrixShell::MatrixShell "callow::MatrixShell::MatrixShell(void *context, const int m, const int
n)

construct with known sizes ";

%feature("docstring")  callow::MatrixShell::~MatrixShell "callow::MatrixShell::~MatrixShell() ";

%feature("docstring")  callow::MatrixShell::set_size "void
callow::MatrixShell::set_size(const int m, const int n=0) ";

%feature("docstring")  callow::MatrixShell::assemble "virtual void
callow::MatrixShell::assemble() ";

%feature("docstring")  callow::MatrixShell::display "virtual void
callow::MatrixShell::display() const ";

%feature("docstring")  callow::MatrixShell::multiply "virtual void
callow::MatrixShell::multiply(const Vector &x, Vector &y)=0 ";

%feature("docstring")  callow::MatrixShell::multiply_transpose "virtual void callow::MatrixShell::multiply_transpose(const Vector &x,
Vector &y)=0 ";


// File: classdetran__geometry_1_1Mesh.xml
%feature("docstring") detran_geometry::Mesh "

Abstract Cartesian mesh class.

Note, the constructors are protected to forbid direct instantiation of
the Mesh class. Rather, use the dimension-specific subclasses. We
could use a pure virtual destructor as an alternative.

C++ includes: Mesh.hh ";

%feature("docstring")  detran_geometry::Mesh::~Mesh "virtual
detran_geometry::Mesh::~Mesh()

Virtual destructor. ";

%feature("docstring")  detran_geometry::Mesh::add_coarse_mesh_map "void detran_geometry::Mesh::add_coarse_mesh_map(std::string map_key,
vec_int mesh_map)

Add map of coarse mesh integer properties.

This is an easy way to set mesh properties for meshes based on simple
coarse mesh regions.

Parameters:
-----------

map_key:  String description of map.

mesh_map:  Logically multi-dimensional map as 1-d vector. ";

%feature("docstring")  detran_geometry::Mesh::add_mesh_map "void
detran_geometry::Mesh::add_mesh_map(std::string map_key, vec_int
mesh_map)

Add map of fine mesh integer properties.

This adds properties for fine meshes directly, and so is meanty for
use with higher level mesh construction, e.g. pin cells, where
assignment is not possible by simple coarse mesh bounds.

If the key exists, this function overwrites the map.

Parameters:
-----------

map_key:  String description of map.

mesh_map:  Logically multi-dimensional map as 1-d vector. ";

%feature("docstring")  detran_geometry::Mesh::number_cells "Mesh::size_t detran_geometry::Mesh::number_cells() const

Return total number of cells. ";

%feature("docstring")  detran_geometry::Mesh::number_cells "Mesh::size_t detran_geometry::Mesh::number_cells(size_t dim) const

Return number of cells in specified dimension. ";

%feature("docstring")  detran_geometry::Mesh::number_cells_x "Mesh::size_t detran_geometry::Mesh::number_cells_x() const

Return number of cells along x axis. ";

%feature("docstring")  detran_geometry::Mesh::number_cells_y "Mesh::size_t detran_geometry::Mesh::number_cells_y() const

Return number of cells along y axis. ";

%feature("docstring")  detran_geometry::Mesh::number_cells_z "Mesh::size_t detran_geometry::Mesh::number_cells_z() const

Return number of cells along z axis. ";

%feature("docstring")  detran_geometry::Mesh::width "double
detran_geometry::Mesh::width(size_t dim, size_t ijk) const

Get cell width in a specified dimension. ";

%feature("docstring")  detran_geometry::Mesh::dx "double
detran_geometry::Mesh::dx(size_t i) const

Get cell width along x axis. ";

%feature("docstring")  detran_geometry::Mesh::dy "double
detran_geometry::Mesh::dy(size_t j) const

Get cell width along y axis. ";

%feature("docstring")  detran_geometry::Mesh::dz "double
detran_geometry::Mesh::dz(size_t k) const

Get cell width along z axis. ";

%feature("docstring")  detran_geometry::Mesh::dx "const Mesh::vec_dbl
& detran_geometry::Mesh::dx() const

Get vector of widths along x axis. ";

%feature("docstring")  detran_geometry::Mesh::dy "const Mesh::vec_dbl
& detran_geometry::Mesh::dy() const

Get vector of widths along y axis. ";

%feature("docstring")  detran_geometry::Mesh::dz "const Mesh::vec_dbl
& detran_geometry::Mesh::dz() const

Get vector of widths along z axis. ";

%feature("docstring")  detran_geometry::Mesh::volume "double
detran_geometry::Mesh::volume(size_t cell) const

Get the cell volume. ";

%feature("docstring")  detran_geometry::Mesh::total_width_x "double
detran_geometry::Mesh::total_width_x() const

Get domain width along x axis. ";

%feature("docstring")  detran_geometry::Mesh::total_width_y "double
detran_geometry::Mesh::total_width_y() const

Get domain width along y axis. ";

%feature("docstring")  detran_geometry::Mesh::total_width_z "double
detran_geometry::Mesh::total_width_z() const

Get domain width along z axis. ";

%feature("docstring")  detran_geometry::Mesh::dimension "Mesh::size_t
detran_geometry::Mesh::dimension() const

Get the mesh dimension. ";

%feature("docstring")  detran_geometry::Mesh::index "Mesh::size_t
detran_geometry::Mesh::index(size_t i, size_t j=0, size_t k=0)

Returns the cardinal index for i, j, and k.

Parameters:
-----------

i:  Index along x axis.

j:  Index along y axis.

k:  Index along z axis.

Cardinal index. ";

%feature("docstring")  detran_geometry::Mesh::cell_to_i "Mesh::size_t
detran_geometry::Mesh::cell_to_i(size_t cell) const

Returns the x index given cardinal index.

Parameters:
-----------

cell:  Cardinal index.

Index along x axis. ";

%feature("docstring")  detran_geometry::Mesh::cell_to_j "Mesh::size_t
detran_geometry::Mesh::cell_to_j(size_t cell) const

Returns the y index given cardinal index.

Parameters:
-----------

cell:  Cardinal index.

Index along y axis. ";

%feature("docstring")  detran_geometry::Mesh::cell_to_k "Mesh::size_t
detran_geometry::Mesh::cell_to_k(size_t cell) const

Returns the z index given cardinal index.

Parameters:
-----------

cell:  Cardinal index.

Index along z axis. ";

%feature("docstring")  detran_geometry::Mesh::mesh_map_exists "bool
detran_geometry::Mesh::mesh_map_exists(std::string map_key)

Check if fine mesh map exists. ";

%feature("docstring")  detran_geometry::Mesh::mesh_map "const
Mesh::vec_int & detran_geometry::Mesh::mesh_map(std::string map_key)

Get map of fine mesh integer properties.

This adds properties for fine meshes directly, and so is meant for use
with higher level mesh construction, e.g. pin cells, where assignment
is not possible by simple coarse mesh bounds.

Parameters:
-----------

m:  Logically multi-dimensional map as 1-d vector. ";

%feature("docstring")  detran_geometry::Mesh::get_mesh_map "const
Mesh::mesh_map_type & detran_geometry::Mesh::get_mesh_map() const

Return a const reference to the full map (useful for IO) ";

%feature("docstring")  detran_geometry::Mesh::display "void
detran_geometry::Mesh::display() const

Display some key features. ";


// File: classdetran__geometry_1_1Mesh1D.xml
%feature("docstring") detran_geometry::Mesh1D "

One-dimensional Cartesian mesh.

This is mostly a convenience interface.

C++ includes: Mesh1D.hh ";

%feature("docstring")  detran_geometry::Mesh1D::Mesh1D "detran_geometry::Mesh1D::Mesh1D(vec_int xfm, vec_dbl xcme, vec_int
mat_map)

Constructor.

Parameters:
-----------

xfm:  Fine meshes per coarse mesh in x dimension.

xcme:  Coarse mesh edges x dimension.

mat_map:  Coarse mesh material map. ";

%feature("docstring")  detran_geometry::Mesh1D::Mesh1D "detran_geometry::Mesh1D::Mesh1D(vec_dbl xfme, vec_int mat_map)

Constructor.

Parameters:
-----------

xfme:  Fine mesh edges x dimension.

mat_map:  Fine mesh material map. ";


// File: classdetran__geometry_1_1Mesh2D.xml
%feature("docstring") detran_geometry::Mesh2D "

Two-dimensional Cartesian mesh.

This is mostly a convenience interface.

C++ includes: Mesh2D.hh ";

%feature("docstring")  detran_geometry::Mesh2D::Mesh2D "detran_geometry::Mesh2D::Mesh2D(vec_int xfm, vec_int yfm, vec_dbl
xcme, vec_dbl ycme, vec_int mat_map)

Constructor.

Parameters:
-----------

xfm:  Fine meshes per coarse mesh in x dimension.

yfm:  Fine meshes per coarse mesh in y dimension.

xcme:  Coarse mesh edges x dimension.

ycme:  Coarse mesh edges y dimension.

mat_map:  Coarse mesh material map. ";

%feature("docstring")  detran_geometry::Mesh2D::Mesh2D "detran_geometry::Mesh2D::Mesh2D(vec_dbl xfme, vec_dbl yfme, vec_int
mat_map)

Constructor.

Parameters:
-----------

xfme:  Fine mesh edges x dimension.

yfme:  Fine mesh edges y dimension.

mat_map:  Fine mesh material map. ";


// File: classdetran__geometry_1_1Mesh3D.xml
%feature("docstring") detran_geometry::Mesh3D "

Three-dimensional Cartesian mesh.

This is mostly a convenience interface.

C++ includes: Mesh3D.hh ";

%feature("docstring")  detran_geometry::Mesh3D::Mesh3D "detran_geometry::Mesh3D::Mesh3D(vec_int xfm, vec_int yfm, vec_int zfm,
vec_dbl xcme, vec_dbl ycme, vec_dbl zcme, vec_int mat_map)

Constructor.

Parameters:
-----------

xfm:  Fine meshes per coarse mesh in x dimension.

yfm:  Fine meshes per coarse mesh in y dimension.

zfm:  Fine meshes per coarse mesh in z dimension.

xcme:  Coarse mesh edges x dimension.

ycme:  Coarse mesh edges y dimension.

zcme:  Coarse mesh edges z dimension.

mat_map:  Coarse mesh material map. ";

%feature("docstring")  detran_geometry::Mesh3D::Mesh3D "detran_geometry::Mesh3D::Mesh3D(vec_dbl xfme, vec_dbl yfme, vec_dbl
zfme, vec_int mat_map)

Constructor.

Parameters:
-----------

xfme:  Fine mesh edges x dimension.

yfme:  Fine mesh edges y dimension.

zfme:  Fine mesh edges z dimension.

mat_map:  Fine mesh material map. ";


// File: classdetran__geometry_1_1MeshMOC.xml
%feature("docstring") detran_geometry::MeshMOC "

Mesh with tracking information.

Todo This is bloat. There should be one mesh object, or perhaps one
geometry object that contains a mesh and tracking

C++ includes: MeshMOC.hh ";

%feature("docstring")  detran_geometry::MeshMOC::MeshMOC "detran_geometry::MeshMOC::MeshMOC(SP_base mesh, SP_trackdb tracks) ";

%feature("docstring")  detran_geometry::MeshMOC::tracks "SP_trackdb
detran_geometry::MeshMOC::tracks() const ";


// File: classdetran_1_1MGCMTSA.xml
%feature("docstring") detran::MGCMTSA "C++ includes: MGCMTSA.hh ";

%feature("docstring")  detran::MGCMTSA::MGCMTSA "detran::MGCMTSA::MGCMTSA(SP_input input, SP_material material, SP_mesh
mesh, SP_scattersource source, size_t cutoff, bool include_fission)

Constructor.

Assuming the within-group transport problem is set up, a KSP object
exists from which the PC is extracted. This PC is passed here to be
constructed and for its application operator to be assigned.

Parameters:
-----------

input:  Input database

material:  Material database

mesh:  Cartesian mesh

source:  Scattering source

cutoff:  First group included in solve ";

%feature("docstring")  detran::MGCMTSA::~MGCMTSA "virtual
detran::MGCMTSA::~MGCMTSA()

virtual destructor ";

%feature("docstring")  detran::MGCMTSA::apply "void
detran::MGCMTSA::apply(Vector &b, Vector &x)

solve Px = b ";

%feature("docstring")  detran::MGCMTSA::multiply "void
detran::MGCMTSA::multiply(const Vector &x, Vector &y) ";

%feature("docstring")  detran::MGCMTSA::multiply_transpose "void
detran::MGCMTSA::multiply_transpose(const Vector &x, Vector &y) ";


// File: classdetran_1_1MGDiffusionSolver.xml
%feature("docstring") detran::MGDiffusionSolver "

Solve a fixed source problem using diffusion.

Fixed source problems can be solved in the following modalities:
fixed, no multiplication

fixed, multiplying via implicit fission (i.e. in the loss matrix)

fixed, multiplying via fission iteration

Mathematically, consider the diffusion equation in operator form:
\\\\[ T \\\\phi = \\\\frac{1}{k} F \\\\phi + Q \\\\, . \\\\] The first
case assumes $ F = 0 $. The second case uses the modified operator $
T' = T - \\\\frac{1}{k} F \\\\phi $. The third case performs the
iteration \\\\[ T \\\\phi^{n+1} = \\\\frac{1}{k}F\\\\phi^{n} + Q \\\\,
. \\\\] This third case is less efficient than the second case, but it
allows one to pull out the solution following each fission iteration.
The number of such fission iterations can be limited by the user.

These three cases are selected via diffusion_fixed_type 0,1,2

C++ includes: MGDiffusionSolver.hh ";

%feature("docstring")  detran::MGDiffusionSolver::MGDiffusionSolver "home robertsj Research detran source src solvers mg MGDiffusionSolver
cc detran::MGDiffusionSolver< D >::MGDiffusionSolver(SP_state state,
SP_material material, SP_boundary boundary, const vec_externalsource
&q_e, SP_fissionsource q_f, bool multiply)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

material:  Material definitions.

boundary:  Boundary fluxes.

external_source:  User-defined external source.

fission_source:  Fission source.

multiply:  Flag for multiplying fixed source problem ";

%feature("docstring")  detran::MGDiffusionSolver::refresh "void
detran::MGDiffusionSolver< D >::refresh()

Refresh the solver. ";

%feature("docstring")  detran::MGDiffusionSolver::lossoperator "SP_lossoperator detran::MGDiffusionSolver< D >::lossoperator()

Return the lossoperator. ";

%feature("docstring")  detran::MGDiffusionSolver::solve "void
detran::MGDiffusionSolver< D >::solve(const double keff=1.0)

Solve the fixed source diffusion problem. ";


// File: classdetran_1_1MGDSA.xml
%feature("docstring") detran::MGDSA "C++ includes: MGDSA.hh ";

%feature("docstring")  detran::MGDSA::MGDSA "home robertsj Research
detran source src solvers mg MGDSA cc detran::MGDSA::MGDSA(SP_input
input, SP_material material, SP_mesh mesh, SP_scattersource source,
size_t cutoff, bool include_fission)

Constructor.

Assuming the within-group transport problem is set up, a KSP object
exists from which the PC is extracted. This PC is passed here to be
constructed and for its application operator to be assigned.

Parameters:
-----------

input:  Input database

material:  Material database

mesh:  Cartesian mesh

source:  Scattering source

cutoff:  First group included in solve ";

%feature("docstring")  detran::MGDSA::~MGDSA "virtual
detran::MGDSA::~MGDSA()

virtual destructor ";

%feature("docstring")  detran::MGDSA::apply "void
detran::MGDSA::apply(Vector &b, Vector &x)

solve Px = b ";

%feature("docstring")  detran::MGDSA::multiply "void
detran::MGDSA::multiply(const Vector &x, Vector &y) ";

%feature("docstring")  detran::MGDSA::multiply_transpose "void
detran::MGDSA::multiply_transpose(const Vector &x, Vector &y) ";


// File: classdetran_1_1MGPreconditioner.xml
%feature("docstring") detran::MGPreconditioner "

Base preconditioner class for multi-group equations.

The diffusion operator only operates on the scalar flux, and so
addition of higher order moments will require restriction and
projection operations.

C++ includes: MGPreconditioner.hh ";

%feature("docstring")  detran::MGPreconditioner::MGPreconditioner "detran::MGPreconditioner::MGPreconditioner(SP_input input, SP_material
material, SP_mesh mesh, size_t cutoff, std::string name=\"MG-PC\")

Constructor.

Parameters:
-----------

input:  Input database

material:  Material database

mesh:  Cartesian mesh

cutoff:  Lowest group to include in the operator ";

%feature("docstring")  detran::MGPreconditioner::~MGPreconditioner "virtual detran::MGPreconditioner::~MGPreconditioner()

virtual destructor ";

%feature("docstring")  detran::MGPreconditioner::apply "virtual void
detran::MGPreconditioner::apply(Vector &b, Vector &x)=0

Solve Px = b. ";


// File: classdetran_1_1MGSolver.xml
%feature("docstring") detran::MGSolver "

Base class for multigroup solvers.

C++ includes: MGSolver.hh ";

%feature("docstring")  detran::MGSolver::MGSolver "detran::MGSolver<
D >::MGSolver(SP_state state, SP_material material, SP_boundary
boundary, const vec_externalsource &q_e, SP_fissionsource q_f, bool
multiply=false)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

material:  Material definitions.

boundary:  Boundary fluxes.

q_e:  Vector of user-defined external sources

q_f:  Fission source.

multiply:  Flag for a multiplying fixed source problem ";

%feature("docstring")  detran::MGSolver::~MGSolver "virtual
detran::MGSolver< D >::~MGSolver()

Virtual destructor. ";

%feature("docstring")  detran::MGSolver::solve "virtual void
detran::MGSolver< D >::solve(const double keff=1.0)=0

Solve the multigroup equations. ";


// File: classdetran_1_1MGSolverGMRES.xml
%feature("docstring") detran::MGSolverGMRES "

Solves the multigroup transport equation via GMRES.

Traditionally, the Gauss-Seidel method has been used for multigroup
problems. For each group, the within-group equation is solved, and the
the fluxes are updated for use in the next group. However, for
problems with significant upscatter, Gauss-Seidel can be quite
expensive, even when GMRES (or some better-than-source-iteration
scheme) is used for the within group solve. As an alternative, we can
apply GMRES (or other Krylov solvers) to the multigroup problem
directly. The linear system is then \\\\[ \\\\left( ( \\\\mathbf{I} -
\\\\left(\\\\begin{array}{ccc} T_1 & \\\\cdots & 0 \\\\\\\\ 0 &
\\\\ddots & 0 \\\\\\\\ 0 & 0 & T_G \\\\end{array}\\\\right) \\\\cdot
\\\\left(\\\\begin{array}{ccc} M & \\\\cdots & 0 \\\\\\\\ 0 &
\\\\ddots & 0 \\\\\\\\ 0 & 0 & M \\\\end{array}\\\\right) \\\\cdot
\\\\left(\\\\begin{array}{ccc} \\\\mathbf{S}_{11} & \\\\cdots &
\\\\mathbf{S}_{1G} \\\\\\\\ \\\\vdots & \\\\ddots & \\\\vdots \\\\\\\\
\\\\mathbf{S}_{G1} & 0 & \\\\mathbf{S}_{GG} \\\\end{array}\\\\right)
\\\\right ) \\\\cdot \\\\left[ \\\\begin{array}{c} \\\\phi_1 \\\\\\\\
\\\\vdots \\\\\\\\ \\\\phi_G \\\\end{array} \\\\right] = \\\\left[
\\\\begin{array}{c} \\\\mathbf{T}_1 q_1 \\\\\\\\ \\\\vdots \\\\\\\\
\\\\mathbf{T}_G q_G \\\\end{array} \\\\right] \\\\, . \\\\] Of course,
this can be written succinctly in the same way as the within- group
equation: \\\\[ (\\\\mathbf{I}-\\\\mathbf{TMS})\\\\phi =
\\\\mathbf{T}q \\\\, , \\\\] where $ \\\\mathbf{T} =
D\\\\mathbf{L}^{-1} $ is the sweeping operator with moment
contributions added implicitly, and where the Krylov vectors are
energy-dependent.

By default, only the energy block in which upscatter occurs is solved
via Krylov methods. Because Gauss-Seidel is exact for downscatter, it
is used for the downscatter-only block. The user can switch this using
\"outer_upscatter_cutoff\".

Reference: Evans, T., Davidson, G. and Mosher, S. \"Parallel
Algorithms for    Fixed-Source and Eigenvalue Problems\", NSTD Seminar
(ORNL), May 27, 2010.

Todo Consider better ways to handle memory between Vec and vector. May
want to devise a moment container based on pointer that allows one to
swap memory temporarily (as done with Vec)

C++ includes: MGSolverGMRES.hh ";

%feature("docstring")  detran::MGSolverGMRES::MGSolverGMRES "home
robertsj Research detran source src solvers mg MGSolverGMRES cc
detran::MGSolverGMRES< D >::MGSolverGMRES(SP_state state, SP_material
material, SP_boundary boundary, const vec_externalsource &q_e,
SP_fissionsource q_f, bool multiply=false)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

material:  Material definitions.

boundary:  Boundary fluxes.

q_e:  Vector of user-defined external sources

q_f:  Fission source.

multiply:  Flag for a multiplying fixed source problem ";

%feature("docstring")  detran::MGSolverGMRES::solve "void
detran::MGSolverGMRES< D >::solve(const double keff=1.0)

Solve the multigroup equations. ";


// File: classdetran_1_1MGSolverGS.xml
%feature("docstring") detran::MGSolverGS "C++ includes: MGSolverGS.hh
";

%feature("docstring")  detran::MGSolverGS::MGSolverGS "home robertsj
Research detran source src solvers mg MGSolverGS cc
detran::MGSolverGS< D >::MGSolverGS(SP_state state, SP_material
material, SP_boundary boundary, const vec_externalsource &q_e,
SP_fissionsource q_f, bool multiply=false)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

material:  Material definitions.

boundary:  Boundary fluxes.

q_e:  Vector of user-defined external sources

q_f:  Fission source.

multiply:  Flag for a multiplying fixed source problem ";

%feature("docstring")  detran::MGSolverGS::solve "void
detran::MGSolverGS< D >::solve(const double keff=1.0)

Solve the multigroup equations. ";


// File: classdetran_1_1MGTransportOperator.xml
%feature("docstring") detran::MGTransportOperator "

Multigroup transport operator.

The multigroup transport operator is defined as ... finish.

C++ includes: MGTransportOperator.hh ";

%feature("docstring")
detran::MGTransportOperator::MGTransportOperator "detran::MGTransportOperator< D >::MGTransportOperator(SP_state state,
SP_boundary boundary, SP_sweeper sweeper, SP_sweepsource source,
size_t cutoff=0)

Constructor.

Parameters:
-----------

state:  state vector

sweeper:  transport sweeper

source:  sweep source

cutoff:  lowest group included in operator ";

%feature("docstring")
detran::MGTransportOperator::~MGTransportOperator "virtual
detran::MGTransportOperator< D >::~MGTransportOperator() ";

%feature("docstring")  detran::MGTransportOperator::display "void
detran::MGTransportOperator< D >::display() const ";

%feature("docstring")  detran::MGTransportOperator::moments_size "size_t detran::MGTransportOperator< D >::moments_size() const ";

%feature("docstring")  detran::MGTransportOperator::boundary_size "size_t detran::MGTransportOperator< D >::boundary_size() const ";

%feature("docstring")  detran::MGTransportOperator::multiply "void
detran::MGTransportOperator< D >::multiply(const Vector &x, Vector &y)
";

%feature("docstring")  detran::MGTransportOperator::multiply_transpose
"virtual void detran::MGTransportOperator< D
>::multiply_transpose(const Vector &x, Vector &y) ";


// File: classdetran_1_1MGTransportSolver.xml
%feature("docstring") detran::MGTransportSolver "

Base class for multigroup transport solvers.

C++ includes: MGTransportSolver.hh ";

%feature("docstring")  detran::MGTransportSolver::MGTransportSolver "detran::MGTransportSolver< D >::MGTransportSolver(SP_state state,
SP_material material, SP_boundary boundary, const vec_externalsource
&q_e, SP_fissionsource q_f, bool multiply)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

material:  Material definitions.

boundary:  Boundary fluxes.

external_source:  User-defined external source.

fission_source:  Fission source.

multiply:  Flag for multiplying fixed source problem ";

%feature("docstring")  detran::MGTransportSolver::~MGTransportSolver "virtual detran::MGTransportSolver< D >::~MGTransportSolver()

Virtual destructor. ";

%feature("docstring")  detran::MGTransportSolver::solve "virtual void
detran::MGTransportSolver< D >::solve(const double keff=1.0)=0

Solve the multigroup equations. ";


// File: classdetran__angle_1_1MomentIndexer.xml
%feature("docstring") detran_angle::MomentIndexer "

Indexex spherical harmonics moments.

The one group (or within-group) scattering source in 3-d is defined
\\\\[ Q(\\\\mathbf{r},\\\\mathbf{\\\\Omega}) = \\\\sum^L_{l=0}
\\\\frac{2l+1}{4\\\\pi} \\\\sum^l_{m=-l} \\\\Sigma_s{l}(\\\\mathbf{r})
\\\\phi^m_l(\\\\mathbf{r}) \\\\: . \\\\] Such sums over $l$ and $m$
are commonplace, in generating sources and also in simply accessing
moments sequentially. This class builds a vector of $(l,m)$ pairs that
can be used in a single loop over all flux moments $\\\\phi^m_l$.

For 3-d, the moments are ordered as \\\\[ [(0,0)] \\\\, , \\\\, \\\\,
\\\\, [(1,-1),\\\\, (1,0)\\\\, (1,1)] \\\\, , \\\\, \\\\, \\\\,
[(2,-2),\\\\,(2,-1)\\\\, \\\\ldots \\\\: . \\\\] For a Legendre order
of $L$, there are $(L+1)^2$ moments.

For 2-d, several moments can be eliminated by symmetry. Only those
moments for which $l+m$ is even are retained (see Hebert). Thus, the
moments are ordered as \\\\[ [(0,0)] \\\\, , \\\\, \\\\, \\\\,
[(1,-1),\\\\, (1,1)] \\\\, , \\\\, \\\\, \\\\, [(2,-2),\\\\,(2,0)\\\\,
\\\\ldots \\\\: . \\\\] For a Legendre order of $L$, there are
$(L+1)(L+2)/2$ moments.

For 1-d, $m=0$, so that the moments are ordered by $l$ alone, from 0
to the given order. Hence, there are $L+1$ moments.

C++ includes: MomentIndexer.hh ";

%feature("docstring")  detran_angle::MomentIndexer::MomentIndexer "home robertsj Research detran source src angle MomentIndexer cc
detran_angle::MomentIndexer::MomentIndexer(const size_t dimension,
const size_t legendre_order)

Constructor.

Parameters:
-----------

dimension:  Problem dimension

legendre_order:  Legendre order of the flux. ";

%feature("docstring")  detran_angle::MomentIndexer::m_index "const
MomentIndexer::vec_int & detran_angle::MomentIndexer::m_index(const
size_t l) const

Return vector of values of m for a given l.

Parameters:
-----------

l:  Legendre degree.

m values ";

%feature("docstring")  detran_angle::MomentIndexer::legendre_order "size_t detran_angle::MomentIndexer::legendre_order() const

Return legendre order. ";

%feature("docstring")  detran_angle::MomentIndexer::number_moments "size_t detran_angle::MomentIndexer::number_moments() const

Return number of moments. ";

%feature("docstring")  detran_angle::MomentIndexer::l "int
detran_angle::MomentIndexer::l(const size_t i) const

Return l value.

Parameters:
-----------

i:  Cardinal moment index.

l value. ";

%feature("docstring")  detran_angle::MomentIndexer::m "int
detran_angle::MomentIndexer::m(const size_t i) const

Return m value.

Parameters:
-----------

i:  Cardinal moment index.

m value. ";

%feature("docstring")  detran_angle::MomentIndexer::index "MomentIndexer::size_t detran_angle::MomentIndexer::index(const size_t
l, const int m) const

Return l value.

Parameters:
-----------

l:  Cardinal moment index.

m:  Cardinal moment index.

Cardinal moment index. ";

%feature("docstring")  detran_angle::MomentIndexer::display "void
detran_angle::MomentIndexer::display() const

Print the indices. ";


// File: classdetran__angle_1_1MomentToDiscrete.xml
%feature("docstring") detran_angle::MomentToDiscrete "

Converts moment-valued unknowns to discrete angle values.

This class defines the operator \\\\[ \\\\mathbf{M} =
\\\\left(\\\\begin{array}{llllllll}
\\\\frac{1}{4\\\\pi}Y^{0}_{0}(\\\\Omega_1) &
\\\\frac{3}{4\\\\pi}Y^{-1}_{1}(\\\\Omega_1) &
\\\\frac{3}{4\\\\pi}Y^{0}_{1}(\\\\Omega_1) &
\\\\frac{3}{4\\\\pi}Y^{1}_{1}(\\\\Omega_1) &
\\\\frac{5}{4\\\\pi}Y^{-2}_{2}(\\\\Omega_1) & \\\\ldots &
\\\\frac{2L+1}{4\\\\pi}Y^{L-1}_{L}(\\\\Omega_1) &
\\\\frac{2L+1}{4\\\\pi}Y^{L}_{L}(\\\\Omega_1) \\\\\\\\
\\\\frac{1}{4\\\\pi}Y^{0}_{0}(\\\\Omega_2) &
\\\\frac{3}{4\\\\pi}Y^{-1}_{1}(\\\\Omega_2) &
\\\\frac{3}{4\\\\pi}Y^{0}_{1}(\\\\Omega_1) &
\\\\frac{3}{4\\\\pi}Y^{1}_{1}(\\\\Omega_2) &
\\\\frac{5}{4\\\\pi}Y^{-2}_{2}(\\\\Omega_2) & \\\\ldots &
\\\\frac{2L+1}{4\\\\pi}Y^{L-1}_{L}(\\\\Omega_2) &
\\\\frac{2L+1}{4\\\\pi}Y^{L}_{L}(\\\\Omega_2) \\\\\\\\ \\\\vdots &
\\\\vdots & \\\\vdots & \\\\vdots & \\\\vdots & & \\\\vdots &
\\\\vdots \\\\\\\\ \\\\frac{1}{4\\\\pi}Y^{0}_{0}(\\\\Omega_{N_n}) &
\\\\frac{3}{4\\\\pi}Y^{-1}_{1}(\\\\Omega_{N_n}) &
\\\\frac{3}{4\\\\pi}Y^{0}_{1}(\\\\Omega_{N_n}) &
\\\\frac{3}{4\\\\pi}Y^{1}_{1}(\\\\Omega_{N_n}) &
\\\\frac{5}{4\\\\pi}Y^{-2}_{2}(\\\\Omega_{N_n}) & \\\\ldots &
\\\\frac{2L+1}{4\\\\pi}Y^{L-1}_{L}(\\\\Omega_{N_n})&
\\\\frac{2L+1}{4\\\\pi}Y^{L}_{L}(\\\\Omega_{N_n}) \\\\\\\\
\\\\end{array}\\\\right) \\\\, . \\\\]

Then, for instance, the angular flux in a particular direction for a
particular cell and group can be approximated by the dot product of
the corresponding row of $\\\\mathbf{M}$ with the vector of flux
moments for that cell and group, i.e. \\\\[ \\\\psi_{i,n} \\\\approx
\\\\frac{1}{4\\\\pi} Y^{0}_{0}(\\\\Omega_n)\\\\phi^{0}_{0} +
\\\\frac{3}{4\\\\pi} Y^{-1}_{1}(\\\\Omega_n)\\\\phi^{-1}_{1} +
\\\\ldots + \\\\frac{2L+1}{4\\\\pi} Y^{L}_{L}(\\\\Omega_n)\\\\phi^L_L
\\\\, , \\\\] for cell $i$ and angle $n$. For a purely isotropic flux,
we have $\\\\phi=\\\\phi^{0}_{0}$ and $\\\\psi_n = \\\\phi/4\\\\pi$ as
we expect (since $ Y^{0}_{0} = 1 $).

For anisotropic scattering of order $L$ in 3-d problems, the number of
columns in $\\\\mathbf{M}$ is $ N_L =(L+1)^2$. In 2-d problems, there
is no variation in the $z$ direction, and hence all values with odd
polar order vanish.

For 1-d problems and scattering order $L$, this reduces to $ N_L =
(L+1) $ and a subsequent simplification of $\\\\mathbf{M}$ into terms
only of Legendre polynomials (and a normalization of $ (2l+1)/2 $
instead of $(2l+1)/4\\\\pi $.

The number of rows $ N_n $ is determined by the number of angles,
which is defined by the quadrature and its order.

The moments of an expanded function are organized as a row of
$\\\\mathbf{M}$, in an array, and one can easily index into that array
for a given $l$ and $m$ using the Moments::index of appropriated
dimension.

See:  Spherical_Harmonics

C++ includes: MomentToDiscrete.hh ";

%feature("docstring")
detran_angle::MomentToDiscrete::MomentToDiscrete "home robertsj
Research detran source src angle MomentToDiscrete cc
detran_angle::MomentToDiscrete::MomentToDiscrete(SP_momentindexer
indexer)

Constructor.

Parameters:
-----------

indexer:  Indexer for spherical harmonic orders ";

%feature("docstring")  detran_angle::MomentToDiscrete::build "void
detran_angle::MomentToDiscrete::build(SP_quadrature q)

Build the moments-to-discrete operator.

Keeping the actual construction outside the constructor allows us to
rebuild the operator for different angular solves using the same
spatial grid. This is useful for coupled forward and adjoint solves,
compact testing of quadrature sets, and potential angular multigrid
schemes.

Parameters:
-----------

q:  Pointer to quadrature ";

%feature("docstring")  detran_angle::MomentToDiscrete::get_row "const
M_Row& detran_angle::MomentToDiscrete::get_row(const size_t angle)
const

Return a row of the operator.

Parameters:
-----------

angle:  Angle i.e. row index.

A row. ";

%feature("docstring")  detran_angle::MomentToDiscrete::row_size "size_t detran_angle::MomentToDiscrete::row_size() const

Return number of moments (length of row in $\\\\mathbf{M}$). ";

%feature("docstring")  detran_angle::MomentToDiscrete::column_size "size_t detran_angle::MomentToDiscrete::column_size() const

Return number of angles (length of column in $\\\\mathbf{M}$). ";

%feature("docstring")  detran_angle::MomentToDiscrete::legendre_order
"size_t detran_angle::MomentToDiscrete::legendre_order() const

Return Legendre expansion order. ";


// File: classdetran_1_1MultiPhysics.xml
%feature("docstring") detran::MultiPhysics "

Container for multiphysics state variables.

This class wraps a map of double vectors that represent various state
variables. The class is \"dumb\" in that it knows nothing about the
problem as is merely a data container.

C++ includes: MultiPhysics.hh ";

%feature("docstring")  detran::MultiPhysics::MultiPhysics "home
robertsj Research detran source src kinetics MultiPhysics cc
detran::MultiPhysics::MultiPhysics(size_t number_variables=0)

Constructor. ";

%feature("docstring")  detran::MultiPhysics::variable "const
MultiPhysics::vec_dbl & detran::MultiPhysics::variable(const size_t
id) const

Const accessor to a physics variable.

Parameters:
-----------

id:  Identifier

Constant reference to physics variable ";

%feature("docstring")  detran::MultiPhysics::variable "MultiPhysics::vec_dbl & detran::MultiPhysics::variable(const size_t
id)

Mutable accessor to a physics variable.

Parameters:
-----------

id:  Identifier

Mutable reference to physics variable ";

%feature("docstring")  detran::MultiPhysics::add_variable "void
detran::MultiPhysics::add_variable(const size_t id, vec_dbl value)

Add a physics variable.

Parameters:
-----------

id:  String identifier

value:  Physics value ";

%feature("docstring")  detran::MultiPhysics::number_variables "size_t
detran::MultiPhysics::number_variables() const

Return the number of physics variables we have. ";

%feature("docstring")  detran::MultiPhysics::display "void
detran::MultiPhysics::display() const

Pretty display of contents. ";


// File: classdetran__utilities_1_1Object.xml
%feature("docstring") detran_utilities::Object "

Abstract class for all objects following DBC.

This class and the associated macros are largely based on the nice
tutorial given at:http://eventhelix.com/realtimemantra/object_oriented
/design_by_contrac t.htm

C++ includes: DBC.hh ";

%feature("docstring")  detran_utilities::Object::~Object "detran_utilities::Object::~Object()

Virtual destructor. ";


// File: classdetran__ortho_1_1OrthogonalBasis.xml
%feature("docstring") detran_ortho::OrthogonalBasis "

Base orthogonal basis class.

C++ includes: OrthogonalBasis.hh ";


// File: classdetran_1_1PC__DSA.xml
%feature("docstring") detran::PC_DSA "

Diffusion synthetic acceleration.

The DSA preconditioning process $ {P}^{-1} $ is defined to be \\\\[
(\\\\mathbf{I} - \\\\mathbf{C}^{-1} \\\\mathbf{S}) \\\\, , \\\\] where
$ \\\\mathbf{C} $ is the one group diffusion operator.

Todo Include fission if treated like scatter

C++ includes: PC_DSA.hh ";

%feature("docstring")  detran::PC_DSA::PC_DSA "home robertsj Research
detran source src solvers wg PC_DSA cc detran::PC_DSA::PC_DSA(SP_input
input, SP_material material, SP_mesh mesh, SP_scattersource source)

Constructor.

Assuming the within-group transport problem is set up, a KSP object
exists from which the PC is extracted. This PC is passed here to be
constructed and for its application operator to be assigned.

Parameters:
-----------

input:  Input database

material:  Material database

mesh:  Cartesian mesh

source:  Scattering source ";

%feature("docstring")  detran::PC_DSA::~PC_DSA "virtual
detran::PC_DSA::~PC_DSA()

virtual destructor ";

%feature("docstring")  detran::PC_DSA::apply "void
detran::PC_DSA::apply(Vector &b, Vector &x)

solve Px = b ";


// File: classcallow_1_1PCIdentity.xml
%feature("docstring") callow::PCIdentity "

Implements an indentity preconditioner (i.e. no preconditioning)

This is mostly for testing purposes.

C++ includes: PCIdentity.hh ";

%feature("docstring")  callow::PCIdentity::PCIdentity "callow::PCIdentity::PCIdentity(double factor=1.0)

Constructor. ";

%feature("docstring")  callow::PCIdentity::~PCIdentity "virtual
callow::PCIdentity::~PCIdentity()

Virtual destructor. ";

%feature("docstring")  callow::PCIdentity::apply "void
callow::PCIdentity::apply(Vector &b, Vector &x)

Solve Px = b. ";


// File: classcallow_1_1PCILU0.xml
%feature("docstring") callow::PCILU0 "

Implements the ILU(0) preconditioner.

Following Saad, ILU(0) is defined

C++ includes: PCILU0.hh ";

%feature("docstring")  callow::PCILU0::PCILU0 "home robertsj Research
detran source src callow preconditioner PCILU0 cc
callow::PCILU0::PCILU0(SP_matrix A)

Construct an ILU0 preconditioner for the explicit matrix A. ";

%feature("docstring")  callow::PCILU0::~PCILU0 "virtual
callow::PCILU0::~PCILU0()

Virtual destructor. ";

%feature("docstring")  callow::PCILU0::apply "void
callow::PCILU0::apply(Vector &b, Vector &x)

Solve Px = b. ";


// File: classcallow_1_1PCJacobi.xml
%feature("docstring") callow::PCJacobi "

Applies a Jacobi preconditioner.

The Jacobi preconditioner is defined by the process \\\\[
\\\\mathbf{P}^{-1} = \\\\mathbf{D}^{-1} =
\\\\mathrm{diag}([a^{-1}_{11}, a^{-1}_{22}, \\\\cdots]^T) \\\\, ,
\\\\] where $ a_{ii} $ is the ith diagonal element of $ \\\\mathbf{A}
$.

If zero elements are found on the diagonal, the inverse is set to the
matrix size.

C++ includes: PCJacobi.hh ";

%feature("docstring")  callow::PCJacobi::PCJacobi "home robertsj
Research detran source src callow preconditioner PCJacobi cc
callow::PCJacobi::PCJacobi(SP_matrix A)

Construct a Jacobi preconditioner for the explicit matrix A. ";

%feature("docstring")  callow::PCJacobi::~PCJacobi "virtual
callow::PCJacobi::~PCJacobi()

Virtual destructor. ";

%feature("docstring")  callow::PCJacobi::apply "void
callow::PCJacobi::apply(Vector &b, Vector &x)

Solve Px = y. ";


// File: classcallow_1_1PCShell.xml
%feature("docstring") callow::PCShell "

Applies a shell preconditioner.

A shell preconditioner allows the user to define their own
preconditioning processes that are potentially matrix free.

C++ includes: PCShell.hh ";

%feature("docstring")  callow::PCShell::PCShell "home robertsj
Research detran source src callow preconditioner PCShell cc
callow::PCShell::PCShell(std::string name=\"PCShell\")

Construct a shell preconditioner. ";

%feature("docstring")  callow::PCShell::~PCShell "virtual
callow::PCShell::~PCShell()

Virtual destructor. ";

%feature("docstring")  callow::PCShell::apply "virtual void
callow::PCShell::apply(Vector &b, Vector &x)=0

Solve Px = b. ";


// File: classdetran__utilities_1_1Point.xml
%feature("docstring") detran_utilities::Point "

Represent a point in three-space.

C++ includes: Point.hh ";

%feature("docstring")  detran_utilities::Point::Point "detran_utilities::Point::Point(double xval=0.0, double yval=0.0,
double zval=0.0) ";

%feature("docstring")  detran_utilities::Point::x "double
detran_utilities::Point::x() const ";

%feature("docstring")  detran_utilities::Point::y "double
detran_utilities::Point::y() const ";

%feature("docstring")  detran_utilities::Point::z "double
detran_utilities::Point::z() const ";


// File: classdetran__angle_1_1PolarQuadrature.xml
%feature("docstring") detran_angle::PolarQuadrature "

Polar quadrature definition, nominally for MOC.

C++ includes: PolarQuadrature.hh ";

%feature("docstring")  detran_angle::PolarQuadrature::PolarQuadrature
"detran_angle::PolarQuadrature::PolarQuadrature(size_t number_polar)

Constructor.

Parameters:
-----------

number_polar:  Number of angles in positive half space ";

%feature("docstring")  detran_angle::PolarQuadrature::~PolarQuadrature
"virtual detran_angle::PolarQuadrature::~PolarQuadrature()

Virtual destructor. ";

%feature("docstring")  detran_angle::PolarQuadrature::number_polar "size_t detran_angle::PolarQuadrature::number_polar() const ";

%feature("docstring")  detran_angle::PolarQuadrature::sin_theta "double detran_angle::PolarQuadrature::sin_theta(size_t p) const

Return a polar sine.

Parameters:
-----------

p:  Polar index in octant ";

%feature("docstring")  detran_angle::PolarQuadrature::cos_theta "double detran_angle::PolarQuadrature::cos_theta(size_t p) const

Return a polar cosine.

Parameters:
-----------

p:  Polar index in octant ";

%feature("docstring")  detran_angle::PolarQuadrature::weight "double
detran_angle::PolarQuadrature::weight(int p) const

Return a polar weight.

Parameters:
-----------

p:  Polar index in octant ";


// File: classcallow_1_1PowerIteration.xml
%feature("docstring") callow::PowerIteration "

Solve the eigenvalue problem with the power method.

C++ includes: PowerIteration.hh ";

%feature("docstring")  callow::PowerIteration::PowerIteration "callow::PowerIteration::PowerIteration(const double tol=1e-6, const
int maxit=100) ";

%feature("docstring")  callow::PowerIteration::~PowerIteration "virtual callow::PowerIteration::~PowerIteration() ";


// File: classcallow_1_1Preconditioner.xml
%feature("docstring") callow::Preconditioner "

Defines a preconditioner for linear solves.

Consider the linear system \\\\[ \\\\mathbf{A}x = b \\\\, . \\\\] When
using iterative methods to solve this system, one can often make the
system easier to solve. Suppose we define an operator $ \\\\mathbf{P}
$ such that $ \\\\mathbf{P}^{-1} \\\\mathbf{A} \\\\approx
\\\\mathbf{I} $. If the action of $ \\\\mathbf{P}^{-1} $ is relatively
easy to compute (as compared to inverting {A}), $ \\\\mathbf{P} $ is a
good preconditioner.

We apply a precondition on the left to obtain the modified system
\\\\[ \\\\mathbf{P^{-1} A}x = \\\\mathbf{P}^{-1} b \\\\, , \\\\] or on
the right to get \\\\[ \\\\mathbf{AP^{-1}}
\\\\overbrace{\\\\mathbf{P}x}^{y} = b \\\\, , \\\\] following which we
solve \\\\[ x = \\\\mathbf{P}^{-1} y \\\\, . \\\\]

Within callow, the Jacobi and ILU(0) preconditioners are available
along with user-defined shell preconditioners. If built with PETSc,
all preconditioners are available (to PETSc) as shells. Otherwise, the
user can set PETSc preconditioners with PetscSolver parameters.

C++ includes: Preconditioner.hh ";

%feature("docstring")  callow::Preconditioner::Preconditioner "callow::Preconditioner::Preconditioner(std::string name) ";

%feature("docstring")  callow::Preconditioner::~Preconditioner "virtual callow::Preconditioner::~Preconditioner() ";

%feature("docstring")  callow::Preconditioner::name "std::string
callow::Preconditioner::name() const

Return the PC name. ";

%feature("docstring")  callow::Preconditioner::apply "virtual void
callow::Preconditioner::apply(Vector &b, Vector &x)=0

solve Px = b ";


// File: classPreconditionerBase.xml
%feature("docstring") PreconditionerBase "

Base class for defining a preconditioner.

Both the within-group and multigroup transport problems can be put in
the form \\\\[ \\\\overbrace{(\\\\mathbf{I} -
\\\\mathbf{D}\\\\mathbf{L}^{-1}\\\\mathbf{MS})}^{\\\\mathbf{A}}
\\\\phi = \\\\overbrace{\\\\mathbf{D} \\\\mathbf{L}^{-1} Q}^{b} \\\\,
, \\\\]

A good preconditioner $ \\\\mathbf{P} $ for this problem should
satisfy \\\\[ \\\\mathbf{P}^{-1} \\\\mathbf{A} \\\\approx
\\\\mathbf{I} \\\\, , \\\\] and, importantly, the action of ${P}^{-1}
$ should be cheap to compute.

Currently, we provide diffusion-based preconditioning that is in
essence equivalent to diffusion synthetic acceleration. Both a within-
group (PreconditionerWG) and a multigroup (PreconditionerMG) are
implemented using the associated within-group diffusion operator
(OneGroupLossOperator) and multigroup diffusion operator
(LossOperator).

In both cases, the preconditioning process ${P}^{-1} $ is defined to
be \\\\[ (\\\\mathbf{I} - \\\\mathbf{C}^{-1} \\\\mathbf{S}) \\\\, ,
\\\\] where $ \\\\mathbf{C} $ is the diffusion operator.

The diffusion operator only operates on the scalar flux, and so
addition of higher order moments will require restriction and
projection operations.

C++ includes: PreconditionerBase.hh ";


// File: classdetran_1_1PreconditionerBase.xml
%feature("docstring") detran::PreconditionerBase "C++ includes:
PreconditionerBase.hh ";

%feature("docstring")  detran::PreconditionerBase::PreconditionerBase
"detran::PreconditionerBase::PreconditionerBase(SP_input input,
SP_material material, SP_mesh mesh, SP_scattersource scattersource) ";

%feature("docstring")  detran::PreconditionerBase::~PreconditionerBase
"virtual detran::PreconditionerBase::~PreconditionerBase()

Virtual destructor. ";

%feature("docstring")  detran::PreconditionerBase::apply "virtual
PetscErrorCode detran::PreconditionerBase::apply(Mat A, Vec x, Vec
y)=0

Apply the preconditioning process, $ \\\\mathbf{P}^{-1} $. ";


// File: classdetran_1_1PreconditionerMG.xml
%feature("docstring") detran::PreconditionerMG "

Multigroup preconditioner.

The multigroup transport problem can be put in the form \\\\[
\\\\overbrace{(\\\\mathbf{I} -
\\\\mathbf{D}\\\\mathbf{L}^{-1}\\\\mathbf{MS})}^{\\\\mathbf{A}}
\\\\phi = \\\\overbrace{\\\\mathbf{D} \\\\mathbf{L}^{-1} Q}^{b} \\\\,
, \\\\]

A good preconditioner $ \\\\mathbf{P} $ for this problem should
satisfy \\\\[ \\\\mathbf{P}^{-1} \\\\mathbf{A} \\\\approx
\\\\mathbf{I} \\\\, , \\\\] and, importantly, the action of ${P}^{-1}
$ should be cheap to compute.

Currently, we provide diffusion-based preconditioner that includes all
in-scatter.

The preconditioning process ${P}^{-1} $ is defined to be \\\\[
(\\\\mathbf{I} - \\\\mathbf{C}^{-1} \\\\mathbf{S}) \\\\, , \\\\] where
$ \\\\mathbf{C} $ is the diffusion operator.

The diffusion operator only operates on the scalar flux, and so
addition of higher order moments will require restriction and
projection operations.

C++ includes: PreconditionerMG.hh ";

%feature("docstring")  detran::PreconditionerMG::PreconditionerMG "detran::PreconditionerMG::PreconditionerMG(SP_input input, SP_material
material, SP_mesh mesh, SP_scattersource source, int
moments_size_group, int boundary_size_group, int upscatter_cutoff, PC
mg_pc)

Constructor.

Assuming the within-group transport problem is set up, a KSP object
exists from which the PC is extracted. This PC is passed here to be
constructed and for its application operator to be assigned.

Parameters:
-----------

input:  Input database

material:  Material database

mesh:  Cartesian mesh

source:  Scattering source

moments_size_group:  Unknown size for group moment vector

boundary_size_group:  Unknown size for group boundary vector

mg_pc:  PETSc PC object to be constructed ";

%feature("docstring")  detran::PreconditionerMG::~PreconditionerMG "detran::PreconditionerMG::~PreconditionerMG()

Destructor. ";

%feature("docstring")  detran::PreconditionerMG::apply "PetscErrorCode detran::PreconditionerMG::apply(Vec x, Vec y)

Apply the preconditioning process, $ \\\\mathbf{P}^{-1} $. ";


// File: classdetran_1_1Precursors.xml
%feature("docstring") detran::Precursors "

Container for precursor concentrations.

C++ includes: Precursors.hh ";

%feature("docstring")  detran::Precursors::Precursors "home robertsj
Research detran source src kinetics Precursors cc
detran::Precursors::Precursors(const size_t number_precursor_groups,
const size_t number_cells)

Constructor.

Parameters:
-----------

number_precursor_groups:  Number of precursor groups

number_cells:  Number of cells ";

%feature("docstring")  detran::Precursors::C "const
Precursors::vec_dbl & detran::Precursors::C(const size_t i) const

Const accessor to a group precursor field.

Parameters:
-----------

i:  Group of precursor requested.

Constant reference to precursor ";

%feature("docstring")  detran::Precursors::C "Precursors::vec_dbl &
detran::Precursors::C(const size_t i)

Mutable accessor to a group moments field.

Parameters:
-----------

i:  Group of precursor requested.

Mutable reference to group moment vector. ";

%feature("docstring")  detran::Precursors::number_precursor_groups "size_t detran::Precursors::number_precursor_groups() const ";

%feature("docstring")  detran::Precursors::number_cells "size_t
detran::Precursors::number_cells() const ";

%feature("docstring")  detran::Precursors::display "void
detran::Precursors::display() const ";


// File: classdetran__angle_1_1ProductQuadrature.xml
%feature("docstring") detran_angle::ProductQuadrature "

Base class for product quadratures.

For the discrete ordinates and characteristics methods, a collocation
in angle is used, and approximations of the integral \\\\[
\\\\int^{\\\\pi}_0 \\\\sin \\\\phi d\\\\phi \\\\int^{2\\\\pi}_{0}
d\\\\theta \\\\, \\\\] are needed over a discrete set of angular
abscissa $ \\\\hat{\\\\Omega}_n $ where \\\\[ \\\\hat{\\\\Omega}_n =
\\\\hat{i} \\\\sin \\\\phi_n \\\\cos \\\\theta_n + \\\\hat{j} \\\\sin
\\\\phi_n \\\\sin \\\\theta_n + \\\\hat{k} \\\\cos \\\\phi_n \\\\]
with associated weights $ w_n $.

A product quadrature separates the angular variables, defining one
dimensional quadrature formulas for each and combining the two as a
product.

As in all quadratures in Detran, we enforce octant symmetry. Hence,
the azimuthal and polar quadratures are defined only within an octant.

Relevant database parameters: quad_number_polar_octant -- number polar
angles per octant

quad_number_azimuth_octant -- number azimuths per octant

C++ includes: ProductQuadrature.hh ";

%feature("docstring")
detran_angle::ProductQuadrature::ProductQuadrature "detran_angle::ProductQuadrature::ProductQuadrature(const size_t dim,
const size_t na, const size_t np, std::string name)

Constructor.

Parameters:
-----------

dim:  Dimension (2 or 3)

na:  Number of azimuths per octant

np:  Number of polar angles per octant

name:   Quadrature name ";

%feature("docstring")
detran_angle::ProductQuadrature::~ProductQuadrature "detran_angle::ProductQuadrature::~ProductQuadrature()=0

Pure virtual destructor. ";

%feature("docstring")  detran_angle::ProductQuadrature::angle "ProductQuadrature::size_t detran_angle::ProductQuadrature::angle(const
size_t a, const size_t p) const

Cardinal angle within octant given azimuth and polar. ";

%feature("docstring")  detran_angle::ProductQuadrature::azimuth "ProductQuadrature::size_t
detran_angle::ProductQuadrature::azimuth(const size_t angle) const

Azimuth index from cardinal within octant. ";

%feature("docstring")  detran_angle::ProductQuadrature::polar "ProductQuadrature::size_t detran_angle::ProductQuadrature::polar(const
size_t angle) const

Polar index from cardinal within octant. ";

%feature("docstring")  detran_angle::ProductQuadrature::sin_theta "double detran_angle::ProductQuadrature::sin_theta(const size_t p)
const ";

%feature("docstring")  detran_angle::ProductQuadrature::cos_theta "double detran_angle::ProductQuadrature::cos_theta(const size_t p)
const ";

%feature("docstring")  detran_angle::ProductQuadrature::phi "double
detran_angle::ProductQuadrature::phi(const size_t a) const ";

%feature("docstring")  detran_angle::ProductQuadrature::sin_phi "double detran_angle::ProductQuadrature::sin_phi(const size_t a) const
";

%feature("docstring")  detran_angle::ProductQuadrature::cos_phi "double detran_angle::ProductQuadrature::cos_phi(const size_t a) const
";

%feature("docstring")
detran_angle::ProductQuadrature::number_azimuths_octant "ProductQuadrature::size_t
detran_angle::ProductQuadrature::number_azimuths_octant() const ";

%feature("docstring")
detran_angle::ProductQuadrature::number_polar_octant "ProductQuadrature::size_t
detran_angle::ProductQuadrature::number_polar_octant() const ";


// File: classdetran_1_1PulsedExternalSource.xml
%feature("docstring") detran::PulsedExternalSource "

External source with Gaussian dependence in time.

Given an external source, a Gaussian shape factor in time is applied.
The Gaussian shape is defined as \\\\[ f(t) = e^{- \\\\frac{(t -
t_{\\\\text{peak}})^2}{2\\\\sigma^2}} \\\\] where the \\\\[ \\\\sigma
= \\\\frac{FWHM}{2\\\\sqrt{2 \\\\ln{2}}} \\\\approx 2.23482 \\\\sigma
\\\\, . \\\\] The user supplies the time at peak, $ t_{\\\\text{peak}}
$ and the full width at half maximum, $ FWHM $.

This shape should be relatively good for modeling pulsed sources.

C++ includes: PulsedExternalSource.hh ";

%feature("docstring")
detran::PulsedExternalSource::PulsedExternalSource "home robertsj
Research detran source src kinetics PulsedExternalSource cc
detran::PulsedExternalSource::PulsedExternalSource(const size_t
number_groups, SP_mesh mesh, SP_externalsource fixed_source, const
double peak_time, const double fwhm, bool discrete=false)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

fixed_source:  Fixed source representing peak pulse

peak_time:  Time of pulse peak

fwhm:  Full width at half maximum of peak

discrete:  Flag to indicate treatment as moment or discrete ";

%feature("docstring")
detran::PulsedExternalSource::~PulsedExternalSource "virtual
detran::PulsedExternalSource::~PulsedExternalSource()

Virtual destructor. ";

%feature("docstring")  detran::PulsedExternalSource::source "double
detran::PulsedExternalSource::source(const size_t cell, const size_t
group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran::PulsedExternalSource::source "double
detran::PulsedExternalSource::source(const size_t cell, const size_t
group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";

%feature("docstring")  detran::PulsedExternalSource::set_time "void
detran::PulsedExternalSource::set_time(const double time)

Set the time of the source.

Parameters:
-----------

time:  Time at which source is evaluated ";


// File: classdetran_1_1PyExecute.xml
%feature("docstring") detran::PyExecute "

Setup and execute the problem from the Python front end.

C++ includes: PyExecute.hh ";

%feature("docstring")  detran::PyExecute::PyExecute "detran::PyExecute< D >::PyExecute(int argc, char *argv[]) ";

%feature("docstring")  detran::PyExecute::initialize "void
detran::PyExecute< D >::initialize(SP_input input, SP_material
material, SP_mesh mesh) ";

%feature("docstring")  detran::PyExecute::solve "void
detran::PyExecute< D >::solve()

Solve the problem. ";

%feature("docstring")  detran::PyExecute::finalize "void
detran::PyExecute< D >::finalize()

Finalize. Closes PETSc, if running, etc. ";


// File: classdetran_1_1PyTimeDependentMaterial.xml
%feature("docstring") detran::PyTimeDependentMaterial "

Base class for time-dependent materials defined in Python.

C++ includes: PyTimeDependentMaterial.hh ";

%feature("docstring")
detran::PyTimeDependentMaterial::PyTimeDependentMaterial "home
robertsj Research detran source src kinetics PyTimeDependentMaterial
cc detran::PyTimeDependentMaterial::PyTimeDependentMaterial(const
size_t number_materials, const size_t number_energy_groups, const
size_t number_precursor_groups, std::string
name=\"PyTimeDependentMaterial\")

Constructor.

Parameters:
-----------

number_materials:  Number of materials

number_energy_groups:  Number of energy groups

number_precursor_groups:  Number of precursor groups

state:   State vector ";

%feature("docstring")
detran::PyTimeDependentMaterial::~PyTimeDependentMaterial "virtual
detran::PyTimeDependentMaterial::~PyTimeDependentMaterial()

Virtual destructor. ";

%feature("docstring")
detran::PyTimeDependentMaterial::set_update_impl "void
detran::PyTimeDependentMaterial::set_update_impl(callback_ptr f, void
*data)

Set the user defined material update. ";


// File: classdetran__angle_1_1Quadrature.xml
%feature("docstring") detran_angle::Quadrature "

Base quadrature class for transport calculations.

For all Detran quadratures, octant symmetry is assumed. As a result,
only the first octant abscissa and weights are stored. This applies to
specializations for product and MOC-specific quadratures.

The octants are ordered as follows, with

*    indices | mu  | eta | xi  *    -------------------------  * 1: N
| +   |  +  | +    (first octant)  *     N+1:2N | -   |  +  | +
(second octant)  *    2N+1:3N | -   |  -  | +    (third octant)  *
3N+1:4N | +   |  -  | +    (fourth octant)  *    4N+1:5N | +   |  +  |
-    ...  *    5N+1:6N | -   |  +  | -  *    6N+1:7N | -   |  -  | - *
7N+1:8N | +   |  -  | -  *

Note that N is the number of angles per quadrant. Outside of the given
pattern, the angles need only be consistently ordered, i.e.

*    abs(mu(i*N+1)) = abs(mu(j*N+1)) for i,j = 0, 1, 2, 3  *

though decreasing absolute value is suggested. This consistency is
helpful for streamlining quadrature sums over all angles.

Because of this ordering, *only the first octant values are stored*.
All other values are found by multiply by the appropriate octant-
dependent sign. This assumes, of course, symmetric quadratures.

C++ includes: Quadrature.hh ";

%feature("docstring")  detran_angle::Quadrature::Quadrature "detran_angle::Quadrature::Quadrature(const size_t dim, const size_t
number_angles, const std::string name)

Constructor.

Parameters:
-----------

dim:  Spatial dimension

number_angles:  Total number of angles

name:  Descriptive name ";

%feature("docstring")  detran_angle::Quadrature::~Quadrature "detran_angle::Quadrature::~Quadrature()=0

Pure virtual destructor. ";

%feature("docstring")  detran_angle::Quadrature::number_angles "Quadrature::size_t detran_angle::Quadrature::number_angles() const

Return total number of angles. ";

%feature("docstring")  detran_angle::Quadrature::number_octants "Quadrature::size_t detran_angle::Quadrature::number_octants() const

Return total number of octants. ";

%feature("docstring")  detran_angle::Quadrature::number_angles_octant
"Quadrature::size_t detran_angle::Quadrature::number_angles_octant()
const

Return number of angles per octant. ";

%feature("docstring")  detran_angle::Quadrature::index "Quadrature::size_t detran_angle::Quadrature::index(const size_t o,
const size_t a)

Return cardinal angle index.

Parameters:
-----------

o:  Octant index

a:  Angle within octant ";

%feature("docstring")  detran_angle::Quadrature::weights "const
Quadrature::vec_dbl & detran_angle::Quadrature::weights() const

Return const reference to weights. ";

%feature("docstring")  detran_angle::Quadrature::cosines "const
Quadrature::vec_dbl & detran_angle::Quadrature::cosines(const size_t
dir) const

Return const reference to a cosine vector.

Parameters:
-----------

dir:  Direction of cosine ";

%feature("docstring")  detran_angle::Quadrature::weight "double
detran_angle::Quadrature::weight(const size_t a) const

Return single weight.

Parameters:
-----------

a:  Angle within octant ";

%feature("docstring")  detran_angle::Quadrature::mu "double
detran_angle::Quadrature::mu(const size_t o, const size_t a) const

Return single $ \\\\mu $.

Parameters:
-----------

o:  Octant index

a:  Angle within octant ";

%feature("docstring")  detran_angle::Quadrature::eta "double
detran_angle::Quadrature::eta(const size_t o, const size_t a) const

Return single $ \\\\eta $.

Parameters:
-----------

o:  Octant index

a:  Angle within octant ";

%feature("docstring")  detran_angle::Quadrature::xi "double
detran_angle::Quadrature::xi(const size_t o, const size_t a) const

Return single $ \\\\xi $.

Parameters:
-----------

o:  Octant index

a:  Angle within octant ";

%feature("docstring")  detran_angle::Quadrature::incident_octant "const Quadrature::vec_int &
detran_angle::Quadrature::incident_octant(const size_t s)

Get vector of incident octants for a side. ";

%feature("docstring")  detran_angle::Quadrature::outgoing_octant "const Quadrature::vec_int &
detran_angle::Quadrature::outgoing_octant(const size_t s)

Get vector of outgoing octants for a side. ";

%feature("docstring")  detran_angle::Quadrature::valid_index "bool
detran_angle::Quadrature::valid_index(const size_t o, const size_t a)
const

Are the indices valid?

Parameters:
-----------

o:  Octant index

a:  Angle within octant ";

%feature("docstring")  detran_angle::Quadrature::dimension "size_t
detran_angle::Quadrature::dimension() const ";

%feature("docstring")  detran_angle::Quadrature::set_adjoint "void
detran_angle::Quadrature::set_adjoint(const bool v=true)

Set adjoint. This changes the octant multipliers. ";

%feature("docstring")  detran_angle::Quadrature::is_adjoint "bool
detran_angle::Quadrature::is_adjoint() const

Is adjoint? ";

%feature("docstring")  detran_angle::Quadrature::name "std::string
detran_angle::Quadrature::name() const ";

%feature("docstring")  detran_angle::Quadrature::display "void
detran_angle::Quadrature::display() const

Pretty print of the first octant parameters. ";


// File: classdetran__angle_1_1QuadratureFactory.xml
%feature("docstring") detran_angle::QuadratureFactory "

Constructs quadratures.

C++ includes: QuadratureFactory.hh ";

%feature("docstring")  detran_angle::QuadratureFactory::build "home
robertsj Research detran source src angle QuadratureFactory cc home
robertsj Research detran source src angle QuadratureFactory cc home
robertsj Research detran source src angle QuadratureFactory cc home
robertsj Research detran source src angle QuadratureFactory cc void
detran_angle::QuadratureFactory::build(SP_quadrature &q, SP_input
input, int dimension)

Build a quadrature set.

Parameters:
-----------

q:  Smart pointer to quadrature being constructed

input:  Smart point to input.

dimension:   Quadrature dimension ";

%feature("docstring")  detran_angle::QuadratureFactory::help "void
detran_angle::QuadratureFactory::help() const

Print out the available quadratures, etc. ";


// File: classdetran__angle_1_1QuadratureMOC.xml
%feature("docstring") detran_angle::QuadratureMOC "

Quadrature specification for MOC transport.

In general, quadrature for MOC has intrinsic space-angle coupling.
Here, spatial aspect is incorporated by defining the start and end
points of each track on a unit cell (i.e. 1x1 cell). These points can
be scaled by the tracker.

Also, several quantities of anticipated use in response generation are
recorded.

Todo It may be useful to define an interface that would require the
polar and azimith components be constructed separately with an
inherited method to build the product set. This would enforce a well-
quantified ordering for future analysis.

C++ includes: QuadratureMOC.hh ";

%feature("docstring")  detran_angle::QuadratureMOC::QuadratureMOC "home robertsj Research detran source src angle QuadratureMOC cc
detran_angle::QuadratureMOC::QuadratureMOC(size_t dim, size_t
num_azimuths_octant, size_t num_polar, std::string name, std::string
polar)

Constructor.

Parameters:
-----------

dim:  Problem dimension (2 or 3d)

num_azimuths_octant:  Number azimuths per octant

num_polar:  Number of polar angles

name:  Name of this azimuthal quadrature

polar:  Name of the requested polar quadrature ";

%feature("docstring")
detran_angle::QuadratureMOC::number_azimuths_octant "size_t
detran_angle::QuadratureMOC::number_azimuths_octant() const ";

%feature("docstring")
detran_angle::QuadratureMOC::number_polar_octant "size_t
detran_angle::QuadratureMOC::number_polar_octant() const ";

%feature("docstring")  detran_angle::QuadratureMOC::number_tracks "size_t detran_angle::QuadratureMOC::number_tracks(size_t a) ";

%feature("docstring")  detran_angle::QuadratureMOC::number_enter "size_t detran_angle::QuadratureMOC::number_enter(size_t a, size_t xy)
const ";

%feature("docstring")  detran_angle::QuadratureMOC::number_exit "int
detran_angle::QuadratureMOC::number_exit(size_t a, size_t xy) const ";

%feature("docstring")  detran_angle::QuadratureMOC::enter "Point
detran_angle::QuadratureMOC::enter(size_t a, size_t i) const ";

%feature("docstring")  detran_angle::QuadratureMOC::exit "Point
detran_angle::QuadratureMOC::exit(size_t a, size_t i) const ";

%feature("docstring")  detran_angle::QuadratureMOC::angle "size_t
detran_angle::QuadratureMOC::angle(size_t a, size_t p) const

Return angle within octant given azimuth and polar. ";

%feature("docstring")  detran_angle::QuadratureMOC::azimuth "size_t
detran_angle::QuadratureMOC::azimuth(size_t angle) const

Return azimuth index for angle within octant. ";

%feature("docstring")  detran_angle::QuadratureMOC::polar "size_t
detran_angle::QuadratureMOC::polar(size_t angle) const

Return polar index for cardinal index. ";

%feature("docstring")  detran_angle::QuadratureMOC::sin_theta "double
detran_angle::QuadratureMOC::sin_theta(size_t p) const ";

%feature("docstring")  detran_angle::QuadratureMOC::cos_theta "double
detran_angle::QuadratureMOC::cos_theta(size_t p) const ";

%feature("docstring")  detran_angle::QuadratureMOC::phi "double
detran_angle::QuadratureMOC::phi(size_t a) const ";

%feature("docstring")  detran_angle::QuadratureMOC::sin_phi "double
detran_angle::QuadratureMOC::sin_phi(size_t a) const ";

%feature("docstring")  detran_angle::QuadratureMOC::cos_phi "double
detran_angle::QuadratureMOC::cos_phi(size_t a) const ";

%feature("docstring")  detran_angle::QuadratureMOC::spacing "double
detran_angle::QuadratureMOC::spacing(size_t a) const ";

%feature("docstring")  detran_angle::QuadratureMOC::azimuth_weight "double detran_angle::QuadratureMOC::azimuth_weight(size_t a) const ";

%feature("docstring")  detran_angle::QuadratureMOC::display_tracks "void detran_angle::QuadratureMOC::display_tracks() const ";


// File: classdetran__angle_1_1QuadrupleRange.xml
%feature("docstring") detran_angle::QuadrupleRange "

Quadruple Range quadrature class.

For two-dimensional transport, Abu-Shumays developed several
quadrature sets that are in some sense analogs of the DPn quadrature.
The sets he developed are product quadratures, which means that polar
and azimuthal angles are independent. To illustrate, consider a
function $ f(\\\\Omega_x,\\\\Omega_y) $. Using a product quadrature,
the integral is approximated as \\\\[ \\\\int_{\\\\Omega}
f(\\\\Omega_x,\\\\Omega_y) d\\\\Omega \\\\approx
\\\\sum^{4N_{\\\\phi}}_{i=1} \\\\sum^{N_{\\\\theta}}_{j=1} w_i w_j
f(\\\\sin \\\\theta_j \\\\cos \\\\phi_i, \\\\sin \\\\theta_j \\\\sin
\\\\phi_i ) \\\\, , \\\\] If a function $ f $ is a polynomial in
$\\\\Omega_x $ and $\\\\Omega_y $ of degree $ L $, a quadrature set
that can exactly integrate $ f $ is denoted order $ L $.

Here, we implement the quadruple range quadrature, which Abu-Shumays
estimated would be well-suited for ``highly heterogeneous problems
with numerous corner singularities''. The quadrature essentially
provides a good approximation for angular fluxes that are represented
by distinct polynomials in $\\\\Omega_x $ and $\\\\Omega_y $ within
each quadrant. We expect such a quadrature to be well-suited for
response function generation, where highly anisotropic fluxes will
result due to incident flux conditions restricted to single surfaces
of a region.

The polar angle quadrature is Abu-Shumays' double-range quadrature,
labeled $ K^d_L $, which satisfies \\\\[ \\\\int^{\\\\pi}_{0}
d\\\\theta \\\\sin \\\\theta \\\\sin^{j} \\\\! \\\\theta =
\\\\int^{1}_{0} x^j \\\\frac{2xdx}{\\\\sqrt{1-x^2}} \\\\, , \\\\]
where $ x = \\\\sin \\\\theta $. To represent terms for $ j = 0
\\\\ldots L $ requires $ L+1 $ equations, sufficient for determining $
(L+1)/2 $ points and $ (L+1)/2 $ weights. Note, this quadrature is
nothing more that a Gaussian quadrature with weight $
w(x)=2x/\\\\sqrt{1-x^2} $ over the range $ 0 \\\\leq x \\\\leq 1 $.
Note furthermore that $ L = 2N_{\\\\theta}-1 $.

The azimuthal angle quadrature is Abu-Shumays' symmetric quadruple-
range quadrature, labeled $ I^{q0}_{L} $. This quadrature satisfies
\\\\[ \\\\int^{\\\\pi/2}_{0} d\\\\phi \\\\cos^n \\\\! \\\\phi
\\\\sin^m \\\\! \\\\phi \\\\approx \\\\sum^{N_{\\\\phi}}_{i=1} w_i
\\\\cos^n \\\\! \\\\phi_i \\\\sin^m \\\\! \\\\phi_i \\\\, . \\\\]
subject to the symmetry constraints \\\\[ \\\\phi_{N_{\\\\phi}+1-i} =
\\\\frac{\\\\pi}{2} - \\\\phi_i \\\\, , \\\\, \\\\, \\\\, \\\\, \\\\,
i = 1 \\\\ldots \\\\frac{N_{\\\\phi}+1}{2} \\\\, \\\\] and \\\\[
w_{N_{\\\\phi}+1-i} = w_i \\\\, , \\\\, \\\\, \\\\, \\\\, \\\\, i = 1
\\\\ldots \\\\frac{N_{\\\\phi}+1}{2} \\\\, . \\\\] Given these
constraints, the moment equations are not all independent, and it
suffices to consider only integrands of form \\\\[ \\\\sin^{4k} \\\\!
\\\\phi\\\\, , \\\\,\\\\, \\\\sin^{2k+1} \\\\! \\\\phi \\\\, , \\\\,
\\\\, \\\\sin^{4k+1} \\\\! \\\\phi \\\\cos \\\\phi \\\\, ,
\\\\,\\\\,\\\\,\\\\,\\\\, k = 0,1\\\\ldots \\\\, . \\\\]

Suppose we select the number of points $ N_{\\\\phi}$, which requires
$ 2N_{\\\\phi} $ equations. The symmetry constraints account for half
of these equations. The remaining $ N_{\\\\phi} $ degrees of freedom
allow for exact integration of polynomials of degree $ L =
N_{\\\\phi}-1 $. Thus, for a given $ L $, we require twice as many
azimuthal points as polar points. When each set has the same $ L $,
they are said to be compatible, i.e. they have the same order of
accuracy. Moreover, it gives a decisive way to choose the number of
polar angles for a given choice of the azimuthal angle count.

For the quadratures as implemented, we use double-precision values as
generated in Maple and Mathematica in tables from which the full set
is generated. In the future, we could also hard code the full set (via
Maple's codegen feature).

The constructor uses the number of azimuthal angles, represented by
the sn_order parameter.

References:

Abu-Shumays, I.K. Nuclear Science and Engineering 64, 299-316 (1977).

C++ includes: QuadrupleRange.hh ";

%feature("docstring")  detran_angle::QuadrupleRange::QuadrupleRange "detran_angle::QuadrupleRange::QuadrupleRange(const size_t order, const
size_t dim=2)

Constructor.

Parameters:
-----------

order:   Quadrature order. ";


// File: classdetran__postprocess_1_1ReactionRates.xml
%feature("docstring") detran_postprocess::ReactionRates "

Computes various reaction rates based on the state.

Once the problem is solved, several quantities are often required for
analysis. These include global net gains and losses, reaction rates in
edit regions such as pins or assemblies, and so on.

C++ includes: ReactionRates.hh ";

%feature("docstring")
detran_postprocess::ReactionRates::ReactionRates "home robertsj
Research detran source src postprocess ReactionRates cc
detran_postprocess::ReactionRates::ReactionRates(SP_material material,
SP_mesh mesh, SP_state state)

Constructor.

Parameters:
-----------

material:  Pin cell pitch (assumed square)

mesh:  Vector of fuel pin radii (can be zero length)

state:  Region material map (cell-center outward) ";

%feature("docstring")
detran_postprocess::ReactionRates::~ReactionRates "virtual
detran_postprocess::ReactionRates::~ReactionRates()

Virtual destructor. ";

%feature("docstring")  detran_postprocess::ReactionRates::region_power
"detran_utilities::vec_dbl
detran_postprocess::ReactionRates::region_power(std::string key,
double scale=1.0)

Relative power of an edit region.

Note, the region powers are returned as a vector in the same order as
the regions are indexed. That means the user is responsible for
mapping the region powers to the end application.

For the built-in pin and assembly arrays, the indexing is natural,
following the same x then y then z ordering as used in Mesh.

Parameters:
-----------

key:  String identifier for the edit region

scale:  Total power used for normalization ";

%feature("docstring")  detran_postprocess::ReactionRates::edit "detran_utilities::vec_dbl
detran_postprocess::ReactionRates::edit(std::string key, const vec_dbl
&fine_mesh_function, bool mean=false)

Edit mesh function.

Given a function defined on the fine mesh, return the function defined
within an edit regions defined by the key. The user can optionally
average the value (rather than just integrate it)

Parameters:
-----------

key:  String identifier for the edit region

fine_mesh_function:  Fine mesh function

mean:  Compute the mean in the edit region ";


// File: classdetran_1_1Reflective.xml
%feature("docstring") detran::Reflective "

Reflective boundary condition so SN problems.

C++ includes: Reflective.hh ";

%feature("docstring")  detran::Reflective::Reflective "detran::Reflective< D >::Reflective(Boundary_T &boundary, const size_t
side, SP_input input, SP_mesh mesh, SP_quadrature quadrature) ";

%feature("docstring")  detran::Reflective::set "void
detran::Reflective< D >::set(const size_t g)

Set initial and/or fixed boundary condition. Reflective does nothing.
";

%feature("docstring")  detran::Reflective::update "void
detran::Reflective< D >::update(const size_t g)

Update a boundary following a sweep. ";

%feature("docstring")  detran::Reflective::update "void
detran::Reflective< D >::update(const size_t g, const size_t o, const
size_t a)

Update a boundary for a given angle following a sweep. ";


// File: classdetran_1_1ReflectiveMOC.xml
%feature("docstring") detran::ReflectiveMOC "

Reflective boundary condition for MOC.

C++ includes: ReflectiveMOC.hh ";

%feature("docstring")  detran::ReflectiveMOC::ReflectiveMOC "detran::ReflectiveMOC< D >::ReflectiveMOC(Boundary_T &boundary, const
size_t side, SP_input input, SP_mesh mesh, SP_quadrature quadrature)
";

%feature("docstring")  detran::ReflectiveMOC::set "void
detran::ReflectiveMOC< D >::set(const size_t g)

Set initial and/or fixed boundary condition. Reflective does nothing.
";

%feature("docstring")  detran::ReflectiveMOC::update "void
detran::ReflectiveMOC< D >::update(const size_t g)

Update a boundary following a sweep. ";

%feature("docstring")  detran::ReflectiveMOC::update "void
detran::ReflectiveMOC< D >::update(const size_t g, const size_t o,
const size_t a)

Update a boundary for a given angle following a sweep. ";


// File: classdetran_1_1ReflectiveSolver.xml
%feature("docstring") detran::ReflectiveSolver "

Solve the reflective boundary condition problem.

This solves \\\\[ \\\\mathbf{L}\\\\psi = q \\\\, , \\\\] with $
\\\\psi $ subject to reflective conditions. In all but pure
reflection, this should take very few sweeps. For pure vacuum
conditions, just one sweep is required.

C++ includes: ReflectiveSolver.hh ";

%feature("docstring")  detran::ReflectiveSolver::ReflectiveSolver "home robertsj Research detran source src solvers ReflectiveSolver cc
detran::ReflectiveSolver< D >::ReflectiveSolver(SP_state state,
SP_boundary boundary, SP_sweeper sweeper, SP_sweepsource source)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

boundary:  Boundary fluxes.

sweeper:

source:  ";

%feature("docstring")  detran::ReflectiveSolver::solve "void
detran::ReflectiveSolver< D >::solve(SP_vector phi)

Solve the reflection equation given flux moments. ";


// File: classcallow_1_1Richardson.xml
%feature("docstring") callow::Richardson "

Uses (modified) Richardson iteration to solve a system.

Richardson iteration solves a linear system via the process \\\\[
x^{(n+1)} = (\\\\mathbf{I - \\\\omega A})x^{(n)} + \\\\omega b \\\\]
where $ \\\\omega $ is something like a relaxation factor that takes
on values between (roughly) 0 and 2. By default, $ \\\\omega = 1 $.

C++ includes: Richardson.hh ";

%feature("docstring")  callow::Richardson::Richardson "callow::Richardson::Richardson(const double atol, const double rtol,
const int maxit, const double omega=1.0) ";

%feature("docstring")  callow::Richardson::~Richardson "virtual
callow::Richardson::~Richardson() ";

%feature("docstring")  callow::Richardson::set_omega "void
callow::Richardson::set_omega(const double omega) ";


// File: classdetran_1_1ScatterSource.xml
%feature("docstring") detran::ScatterSource "

Methods for constructing various scattering sources.

Todo Implement something in group bounds for adjoint

C++ includes: ScatterSource.hh ";

%feature("docstring")  detran::ScatterSource::ScatterSource "detran::ScatterSource::ScatterSource(SP_mesh mesh, SP_material
material, SP_state state)

Constructor.

Parameters:
-----------

mesh:  Pointer to mesh

material:  Pointer to material

state:  Pointer to state vector ";

%feature("docstring")
detran::ScatterSource::build_within_group_source "void
detran::ScatterSource::build_within_group_source(const size_t g, const
moments_type &phi, moments_type &source)

Build the within group scattering source.

This constructs \\\\[ q_g = \\\\mathbf{S}_{gg} \\\\phi_g \\\\, . \\\\]

Parameters:
-----------

g:  Group for this problem

phi:  Const reference to group flux.

source:  Mutable reference to moments source. ";

%feature("docstring")  detran::ScatterSource::build_in_scatter_source
"void detran::ScatterSource::build_in_scatter_source(const size_t g,
moments_type &source)

Build the in-scatter source.

This constructs \\\\[ q_g = \\\\sum^G_{g',g\\\\ne g'}
\\\\mathbf{S}_{gg'}\\\\phi_{g'} \\\\, . \\\\]

This *assumes* the state is up-to-date.

Parameters:
-----------

g:  Group for this problem

source:  Mutable reference to moments source.

g_down:  Highest group to contribute to downscatter

g_up:  Lowest group to contribute to upscatter ";

%feature("docstring")  detran::ScatterSource::build_downscatter_source
"void detran::ScatterSource::build_downscatter_source(const size_t g,
const size_t g_cutoff, moments_type &source)

Build the downscatter source.

This constructs \\\\[ q_g = \\\\sum^{g_{cutoff}_{g'}
\\\\mathbf{S}_{gg'}\\\\phi_{g'} \\\\, . \\\\]

This assumes the state is up-to-date.

This is useful when creating the fixed source for Krylov multigroup
solves when Gauss-Seidel has been used for the downscatter block.

Parameters:
-----------

g:  Group for this problem

g_down:  Highest group to contribute to downscatter

source:  Mutable reference to moments source. ";

%feature("docstring")  detran::ScatterSource::build_total_group_source
"void detran::ScatterSource::build_total_group_source(const size_t g,
const size_t g_cutoff, const State::vec_moments_type &phi,
moments_type &source)

Build the total scatter source.

In some cases, including all scattering is required, as is the case
when performing multigroup Krylov solves. This constructs

\\\\[ q_g = \\\\sum^G_{g'} \\\\mathbf{S}_{gg'}\\\\phi_{g'} \\\\, .
\\\\]

Because Gauss-Seidell can be used to solve downscatter blocks, a
cutoff group is passed to exclude the solved portion of the problem.

Parameters:
-----------

g:  Group for this problem.

g_cutoff:  Highest group to contribute to downscatter.

phi:  Const reference to multigroup flux moments.

source:  Mutable reference to moments source. ";


// File: classdetran__geometry_1_1Segment.xml
%feature("docstring") detran_geometry::Segment "

One segment along a track.

Each segment traverses one flat source region. To characterize the
segment, we need its length and an index to the region.

Todo For acceleration purposes, we might also need to know what coarse
mesh boundary it intersects, if any.

C++ includes: Segment.hh ";

%feature("docstring")  detran_geometry::Segment::Segment "detran_geometry::Segment::Segment(const size_t r, double l) ";

%feature("docstring")  detran_geometry::Segment::region "int
detran_geometry::Segment::region() const

Return the flat source region. ";

%feature("docstring")  detran_geometry::Segment::length "double
detran_geometry::Segment::length() const

Return the segment length. ";

%feature("docstring")  detran_geometry::Segment::scale "void
detran_geometry::Segment::scale(double v)

Scale a segment. ";


// File: classdetran__ioutils_1_1SiloOutput.xml
%feature("docstring") detran_ioutils::SiloOutput "

Write mesh data to a Silo file.

C++ includes: SiloOutput.hh ";

%feature("docstring")  detran_ioutils::SiloOutput::SiloOutput "home
robertsj Research detran source src ioutils SiloOutput cc
detran_ioutils::SiloOutput::SiloOutput(SP_mesh mesh)

Constructor.

Parameters:
-----------

input:  Input database

mesh:  Cartesian mesh ";

%feature("docstring")  detran_ioutils::SiloOutput::~SiloOutput "detran_ioutils::SiloOutput::~SiloOutput()

Destructor. ";

%feature("docstring")  detran_ioutils::SiloOutput::initialize "bool
detran_ioutils::SiloOutput::initialize(const std::string filename)

Initialize file. ";

%feature("docstring")  detran_ioutils::SiloOutput::write_mesh_map "bool detran_ioutils::SiloOutput::write_mesh_map(const std::string
&key)

Write a mesh map to file.

Parameters:
-----------

key:  Key of the map to write

True for successful write ";

%feature("docstring")  detran_ioutils::SiloOutput::write_scalar_flux "bool detran_ioutils::SiloOutput::write_scalar_flux(SP_state state)

Write the multigroup scalar flux moments to file.

Parameters:
-----------

state:  State vector container

True for successful write ";

%feature("docstring")  detran_ioutils::SiloOutput::write_angular_flux
"bool detran_ioutils::SiloOutput::write_angular_flux(SP_state state,
SP_quadrature quad)

Write the angular flux to file.

Parameters:
-----------

state:  State vector container

True for successful write ";

%feature("docstring")  detran_ioutils::SiloOutput::write_time_flux "bool detran_ioutils::SiloOutput::write_time_flux(const int step,
SP_state state, bool do_psi)

Write the time dependent fluxes to file.

Parameters:
-----------

step:  Time step

state:  State vector container

True for successful write ";

%feature("docstring")  detran_ioutils::SiloOutput::finalize "void
detran_ioutils::SiloOutput::finalize()

Close file. ";

%feature("docstring")  detran_ioutils::SiloOutput::add_material "void
detran_ioutils::SiloOutput::add_material(SP_material material)

Add the material database.

The user can create maps of the cross section values.

Parameters:
-----------

material:  Material database ";


// File: classdetran_1_1Solver.xml
%feature("docstring") detran::Solver "

Boilerplate shared by all solvers.

C++ includes: Solver.hh ";

%feature("docstring")  detran::Solver::Solver "detran::Solver< D
>::Solver(SP_state state, SP_material material, SP_boundary boundary,
const vec_externalsource &q_e, SP_fissionsource q_f)

Constructor.

Parameters:
-----------

input:  Input database.

state:   State vectors, etc.

mesh:  Problem mesh.

material:  Material definitions.

boundary:  Boundary fluxes.

q_e:  Vector of user-defined external sources

q_f:  Fission source.

multiply:  Flag for a multiplying fixed source problem ";

%feature("docstring")  detran::Solver::~Solver "detran::Solver< D
>::~Solver()=0

Pure virtual destructor. ";

%feature("docstring")  detran::Solver::set_tolerance "void
detran::Solver< D >::set_tolerance(double tol)

Reset the tolerance. ";

%feature("docstring")  detran::Solver::set_max_iters "void
detran::Solver< D >::set_max_iters(int max_iters)

Reset the maximum iterations. ";

%feature("docstring")  detran::Solver::refresh "virtual void
detran::Solver< D >::refresh()

Refresh the solver.

This can be reimplemented by a subclass if internal structures must be
updated due to some change (e.g. material perturbation) ";


// File: classdetran__utilities_1_1SP.xml
%feature("docstring") detran_utilities::SP "

Smart pointer implementation that does reference counting.

The smart pointer provides a \"safe\" encapsulation for a standard C++
pointer. Consider: A function new's an object and return the pointer
to that object as its return value. Now it is the caller's
responsibility to free the object. What if the caller passes the
pointer to other objects or functions? What if it is not known which
will be deleted first or last?

Instead the function can return a \"smart pointer\". This SP class
uses reference counting to determine the number of current users of a
pointer. Each time an SP goes out of scope, the reference count is
decremented. When the last user of a pointer is done, the pointer is
freed.

Note: I am calling this a \"smart pointer\", not a \"safe pointer\".
There are clearly ways you can hose this. In particular, when you bind
an SP<T> to a T*, you yield all rights to the T*. You'd better not
squirrel the bare pointer away somewhere and expect to clandestinely
use it in other ways or places--death will be sure to follow.
Consequently then, the safest way to use this smart pointer, is to
bind it to the contained pointer and then always use the smart
pointer. Immediately returning the smart pointer as a return value,
allowing the original bare pointer to go out of scope never to be seen
again, is one good example of how to use this.

One good example of bad usage is assigning the same dumb pointer to
multiple SPs. Consider: Unfortunately, there is no way to check if
another SP owns the dumb pointer that you give to a SP. This is simply
something that needs to be watched by the programmer.

Note, SP is minimally thread-safe in that worker threads can make
copies of the SP (which increments the counter) and when out of scope,
the counter is decremented, all safely in the counter. The objects to
which these SP's point are not thread safe, but there are few, if any,
cases where that behavior would be required.

C++ includes: SP.hh ";

%feature("docstring")  detran_utilities::SP::SP "detran_utilities::SP< T >::SP()

Default constructor. ";

%feature("docstring")  detran_utilities::SP::SP "detran_utilities::SP< T >::SP(T *p_in)

Explicit constructor for type T *.

This constructor is used to initialize a SP with a pointer, ie. Once a
pointer is \"given\" to a SP, the SP takes control. This means that
spf will delete f when the last SP to f is destroyed.

Parameters:
-----------

p_in:  pointer to type T ";

%feature("docstring")  detran_utilities::SP::SP "detran_utilities::SP< T >::SP(X *px_in)

Explicit constructor for type X *.

This constructor is used to initialize a base class smart pointer of
type T with a derived class pointer of type X or, equivalently, any
types in which X * is convertible to T * through a dynamic_cast.
Consider, The pointer to X must not be equal to NULL. The SP owns the
pointer when it is constructed.

Parameters:
-----------

px_in:  pointer to type X that is convertible to T * ";

%feature("docstring")  detran_utilities::SP::SP "detran_utilities::SP< T >::SP(const SP< T > &sp_in)

Copy constructor for SP<T>.

Parameters:
-----------

sp_in:  smart pointer of type SP<T> ";

%feature("docstring")  detran_utilities::SP::SP "detran_utilities::SP< T >::SP(const SP< X > &spx_in)

Copy constructor for SP<X>.

This copy constructor Requires that X * is convertible to T * through
a dynamic_cast. The pointer in spx_in can point to NULL; however, it
must still be convertible to T * through a dynamic_cast.

Parameters:
-----------

spx_in:  smart pointer of type SP<X> ";

%feature("docstring")  detran_utilities::SP::~SP "detran_utilities::SP< T >::~SP()

Destructor, memory is released when count goes to zero. ";

%feature("docstring")  detran_utilities::SP::bp "T*
detran_utilities::SP< T >::bp() const

Get the base-class pointer; better know what you are doing. ";


// File: classdetran__angle_1_1SphericalHarmonics.xml
%feature("docstring") detran_angle::SphericalHarmonics "

Spherical harmonics generation for anisotropic scattering.

In 2-d and 3-d problems, an expansion of the angular flux in spherical
harmonics is required to treat anisotropic scattering. In 1-d
problems, the angular flux is expanded in Legendre polynomials, which
can be considered to be a subset of the spherical harmonics.

We follow the presentation of Hebert.

A function of a directional cosine $\\\\xi$, $ f(\\\\xi) $, can be
represented as in $L$-th order Legendre expansion: \\\\[ f( \\\\xi)
\\\\approx \\\\sum^L_{l=0} \\\\frac{2l+1}{2}P_l(\\\\xi)f_l \\\\, ,
\\\\] where the Legendre coefficients are defined \\\\[ f_l =
\\\\int^{1}_{-1} d\\\\xi P_l(\\\\xi) f(\\\\xi) \\\\, . \\\\] In the
standard 1-d formulation, in which the azimuthal dependence is removed
by integration, the zeroth Legendre moment angular flux $
\\\\psi(z,\\\\xi) $ is \\\\[ \\\\phi_0(z) = \\\\int^{1}_{-1} d\\\\xi
\\\\psi(z,\\\\xi) \\\\, , \\\\] which we recognize as the scalar flux.
The scattering source is defined in 1D as \\\\[ Q(z,\\\\xi) =
\\\\int^{2\\\\pi}_{0} d\\\\phi' \\\\int^{1}_{-1}d\\\\xi'
\\\\Sigma_s(z,\\\\mu_0)\\\\psi(z,\\\\xi') \\\\\\\\ \\\\] and can
likewise be expanded (skipping several important steps! see e.g. the
22106 notes) as \\\\[ Q(z,\\\\xi) \\\\approx \\\\sum^{L}_{l=0}
\\\\frac{2l+1}{2}\\\\Sigma_{sl}(z)P_l(\\\\xi)\\\\phi_l(z) \\\\, ,
\\\\] where the Legendre moments $ \\\\Sigma_{sl} $ are provided in a
cross-section library. Note, for consistency with the multidimensional
coordinates, we place 1D problems along the $z$-axis. For 1D problems
using an $L$-th order expansion (for scattering), there are $L+1$ flux
moments per node (and possibly more than one node per mesh).

For 2D and 3D problems, the azimuthal dependence requires use of the
spherical harmonics for expansions.

The spherical harmonics $Y^m_l(\\\\Omega) = Y^m_l(\\\\xi,\\\\varphi)$
(in the source, denoted $Y_{lm}$ ) are defined \\\\[
Y^m_l(\\\\xi,\\\\varphi) = \\\\sqrt{ (2-\\\\delta_{m,0})
\\\\frac{(l-|m|)!}{(l+|m|)!} } P^{|m|}_l (\\\\xi) \\\\mathcal{T}_m
(\\\\varphi ) \\\\ . \\\\] where the associated Legendre polynomials
are defined in terms of Legendre polynomials as \\\\[ P^m_l(\\\\xi) =
(1-\\\\xi^2)^{m/2} \\\\frac{d^m}{d\\\\xi^m} P_l(\\\\xi) \\\\, ,
\\\\,\\\\,\\\\,\\\\,\\\\, m \\\\geq 0 \\\\, , \\\\] and the
trigonometric functions are \\\\[ \\\\mathcal{T}_m (\\\\varphi )
\\\\left\\\\{ \\\\begin{array}{l l} \\\\cos(m\\\\varphi) & \\\\quad
\\\\text{if $m \\\\geq 0$ }\\\\\\\\ \\\\sin(|m|\\\\varphi) & \\\\quad
\\\\text{otherwise} \\\\, . \\\\\\\\ \\\\end{array} \\\\right. \\\\]
The associated Legendre polynomials are as given use the so-called
Ferrer definition, omitting the typically standard factor of $ (-1)^m
$. Hebert notes this representation is helpful since \\\\[
\\\\mathbf{\\\\Omega} = \\\\begin{pmatrix} \\\\mu \\\\\\\\ \\\\eta
\\\\\\\\ \\\\xi \\\\end{pmatrix} = \\\\begin{pmatrix}
\\\\sin(\\\\theta)\\\\cos(\\\\varphi) \\\\\\\\
\\\\sin(\\\\theta)\\\\sin(\\\\varphi) \\\\\\\\ \\\\cos(\\\\theta)
\\\\end{pmatrix} = \\\\begin{pmatrix} Y^{1}_{1}(\\\\xi,\\\\varphi)
\\\\\\\\ Y^{-1}_{1}(\\\\xi,\\\\varphi) \\\\\\\\
Y^{0}_{1}(\\\\xi,\\\\varphi) \\\\end{pmatrix} \\\\, . \\\\] The
spherical harmonics through $L=3$ in terms of the directional cosines
are \\\\[ Y^{0}_{0} = 1 \\\\, , \\\\] \\\\[ Y^{-1}_{1}=\\\\eta \\\\, ,
\\\\quad Y^{0}_{1}=\\\\xi \\\\, , \\\\quad Y^{1}_{1}=\\\\mu \\\\]
\\\\[ Y^{-2}_{2} = \\\\sqrt{3}\\\\mu\\\\eta \\\\, , \\\\quad
Y^{-1}_{2}=\\\\sqrt{3}\\\\xi\\\\eta \\\\, , \\\\quad
Y^{0}_{2}=\\\\frac{1}{2}(3\\\\xi^2-1) \\\\, , \\\\quad
Y^{1}_{2}=\\\\sqrt{3}\\\\xi\\\\mu \\\\, , \\\\quad
Y^{2}_{2}=\\\\frac{\\\\sqrt{3}}{2}(\\\\mu^2-\\\\eta^2) \\\\, . \\\\]
\\\\[ Y^{-3}_{3} =
\\\\sqrt{\\\\frac{5}{8}}\\\\eta(3\\\\mu^2-\\\\eta^2) \\\\, , \\\\quad
Y^{-2}_{3} = \\\\sqrt{15}\\\\xi\\\\mu\\\\eta \\\\, , \\\\quad
Y^{-1}_{3} = \\\\sqrt{\\\\frac{3}{8}}\\\\eta(5\\\\xi^2-1) \\\\, ,
\\\\quad Y^{0}_{3} = \\\\frac{1}{2}(5\\\\xi^3-3\\\\xi) \\\\, ,
\\\\quad Y^{1}_{3} = \\\\sqrt{\\\\frac{3}{8}}\\\\mu(5\\\\xi^2-1) \\\\,
, \\\\quad Y^{2}_{3} =
\\\\sqrt{\\\\frac{15}{4}}\\\\xi(\\\\mu^2-\\\\eta^2) \\\\, , \\\\quad
Y^{3}_{3} = \\\\sqrt{\\\\frac{5}{8}}\\\\mu(\\\\mu^2-3\\\\eta^2) \\\\,
. \\\\] The trigonometric functions, associated Legendre polynomials,
and spherical harmonics satisfy the following orthogonality
conditions: \\\\[ \\\\int^{\\\\pi}_{-\\\\pi} d\\\\varphi
\\\\mathcal{T}_n (\\\\varphi ) \\\\mathcal{T}_{n'}(\\\\varphi ) =
\\\\pi(1+\\\\delta_{n,0})\\\\delta_{n,n'} \\\\, , \\\\] \\\\[
\\\\int^{1}_{-1} d\\\\mu P^m_l (\\\\mu) P^{m'}_{l'} =
\\\\frac{2(l+m)!}{(2l+1)(l-m)!}\\\\delta_{l,l'} \\\\, , \\\\] and
\\\\[ \\\\int_{4\\\\pi} d^2\\\\Omega
Y^m_l(\\\\mathbf{\\\\Omega})Y^m_{l'}(\\\\mathbf{\\\\Omega}) =
\\\\int^{2\\\\pi}_{0} d\\\\phi' \\\\int^{1}_{-1}d\\\\mu'
Y^m_l(\\\\theta,\\\\varphi) Y^{m'}_{l'}(\\\\theta,\\\\varphi) =
\\\\frac{4\\\\pi}{2l+1} \\\\delta_{l,l'}\\\\delta_{m,m'} \\\\, . \\\\]
Then, angular flux is expanded following \\\\[
\\\\psi(\\\\mathbf{r},\\\\mathbf{\\\\Omega}) \\\\approx
\\\\sum^{L}_{l=0} \\\\frac{2l+1}{4\\\\pi} \\\\sum^{l}_{m=-l} Y^m_l
(\\\\mathbf{\\\\Omega}) \\\\phi^m_l(\\\\mathbf{r}) \\\\, , \\\\] where
\\\\[ \\\\phi^m_l(\\\\mathbf{r}) = \\\\int_{4\\\\pi} d^2\\\\Omega
Y^m_l(\\\\mathbf{\\\\Omega})
\\\\psi(\\\\mathbf{r},\\\\mathbf{\\\\Omega}) \\\\, . \\\\] The
additional $2\\\\pi$ in the expansion comes from the isotropy in the
azimuthal angle, multiplied away in 1D (a 1D angular flux is really in
units of 1/rad, not 1/steradian).

Again, we note that the zeroth order moment is the angular flux. For
an $L$th order expansion, there are $(L+1)^2$ moments.

For 2-D, the moments $\\\\phi^0_{l>0}$ vanish identically. This is
because a 2D problem is defined such that there is no variation in the
polar ( $z$) direction. The moments $\\\\phi^0_{l>0}$ represent net
changes in the $z$ direction, and so should be eliminated for
efficiency. For 2D problems, the number of moments is therefore
$(L+1)^2-L$. This is handled in Moment_to_Discrete.  It should be
noted that the for $m=0 $, the spherical harmonics as defined reduce
to the Legendre polynomials in $\\\\xi$. Hence, this class can be used
for 1D expansions.

Currently, only the spherical harmonics through $ L=3 $ are
implemented, though adding capability for arbitrary orders should be
straightforward if external libraries are used.

Legendre expansions are typically referred to in terms of the
\"order\", which we have denoted via L. In many mathematics texts, the
term \"degree\" is used instead. Here, both \"Legendre order\" and
\"spherical harmonic degree\" will refer to the l subscript, while the
m subscript represents the \"spherical harmonic order\".

Alain Hebert, Applied Reactor Physics, Presses Internationales
Polytechnique, Montreal, 2009.

C++ includes: SphericalHarmonics.hh ";


// File: classdetran__utilities_1_1SPref.xml
%feature("docstring") detran_utilities::SPref "

Reference counter for SP class.

This reference counter is thread safe, and allows SP's to be used as
pointers to unmutable objects.

C++ includes: SP.hh ";

%feature("docstring")  detran_utilities::SPref::SPref "detran_utilities::SPref::SPref(int r=1)

Constructor. ";

%feature("docstring")  detran_utilities::SPref::refs "int
detran_utilities::SPref::refs() const

Return the reference count. ";

%feature("docstring")  detran_utilities::SPref::increment "void
detran_utilities::SPref::increment()

Increment the reference count. ";

%feature("docstring")  detran_utilities::SPref::decrement "void
detran_utilities::SPref::decrement()

Decrement the reference count. ";


// File: classdetran_1_1State.xml
%feature("docstring") detran::State "

Problem state.

The state of the problem at any point in the solution process can be
defined by the unknown flux moments (scalar flux and higher moments),
as these quantities are sufficient to describe reaction rates, which
is typically what we need (e.g. doses or fission rates). For
eigenvalue problems, keff is also important.

When needed, the underlying data array for a group can be accessed as
demonstrated by the following (I think!). This would be useful e.g.
for filling a PETSc Vec object to use in one of their Krylov schemes
(rather than copying into another array).

Todo Test whether referencing a Moments_Field object at [0] actually
yields the array underneath.

Relevant input entries: number_groups (int)

store_angular_flux (int)

C++ includes: State.hh ";

%feature("docstring")  detran::State::State "detran::State::State(SP_input input, SP_mesh mesh, SP_quadrature
quadrature=SP_quadrature(0))

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

quadrature:  Angular quadrature. ";

%feature("docstring")  detran::State::phi "const State::moments_type
& detran::State::phi(const size_t g) const

Const accessor to a group moments field.

Parameters:
-----------

g:  Group of field requested.

Constant reference to group moment vector. ";

%feature("docstring")  detran::State::phi "State::moments_type &
detran::State::phi(const size_t g)

Mutable accessor to a group moments field.

This is to be used for copying, i.e.

Parameters:
-----------

g:  Group of field requested.

Mutable reference to group moment vector. ";

%feature("docstring")  detran::State::all_phi "const
State::group_moments_type & detran::State::all_phi() const

Const access to all group moments. ";

%feature("docstring")  detran::State::all_phi "State::group_moments_type & detran::State::all_phi()

Mutable access to all group moments. ";

%feature("docstring")  detran::State::set_moments "void
detran::State::set_moments(const size_t g, std::vector< double > &f)

Set a group moment vector.

Parameters:
-----------

g:  Energy group

f:  User-defined moment vector. ";

%feature("docstring")  detran::State::psi "const angular_flux_type&
detran::State::psi(const size_t g, const size_t o, const size_t a)
const

Const accessor to a group angular flux.

Parameters:
-----------

g:  Group of field requested.

o:  Octant

a:  Angle within octant

Constant reference to group angular flux vector. ";

%feature("docstring")  detran::State::psi "State::angular_flux_type &
detran::State::psi(const size_t g, const size_t o, const size_t a)

Mutable accessor to a group angular flux.

Parameters:
-----------

g:  Group of field requested.

o:  Octant

a:  Angle within octant

Mutable reference to group angular flux vector. ";

%feature("docstring")  detran::State::eigenvalue "double
detran::State::eigenvalue() const ";

%feature("docstring")  detran::State::set_eigenvalue "void
detran::State::set_eigenvalue(const double v) ";

%feature("docstring")  detran::State::get_input "SP_input
detran::State::get_input() ";

%feature("docstring")  detran::State::get_mesh "SP_mesh
detran::State::get_mesh() ";

%feature("docstring")  detran::State::get_momentindexer "SP_momentindexer detran::State::get_momentindexer() ";

%feature("docstring")  detran::State::moments_size "size_t
detran::State::moments_size() const ";

%feature("docstring")  detran::State::get_quadrature "SP_quadrature
detran::State::get_quadrature() ";

%feature("docstring")  detran::State::number_groups "size_t
detran::State::number_groups() const ";

%feature("docstring")  detran::State::store_angular_flux "bool
detran::State::store_angular_flux() const ";

%feature("docstring")  detran::State::clear "void
detran::State::clear()

Zero out the state. ";

%feature("docstring")  detran::State::scale "void
detran::State::scale(const double f)

Scale the state by a constant. ";

%feature("docstring")  detran::State::display "void
detran::State::display() const

Format display of flux. ";


// File: classdetran_1_1StupidParser.xml
%feature("docstring") detran::StupidParser "

Parse a really simple input file.

Because I'm exceedingly lazy and want to use Python for all my
\"input\" now that I've tried it, I really want to spend little time
on a text-based input format. Hence, I'm doing the absolutely easiest
(IMHO) thing possible. For those who need more flexibility, write a
wrapper or better, contribute a more reasonable parser. The input file
will be a series of lines similar to the following, which is a simple
working example

# General input (the \"#\" is ignored)    int      number_groups 2
int      dimension               2    dbl inner_max_tolerance     1e-4
str      problem_type            fixed # Mesh (blank lines also
ignored)    mesh mesh_xcme 0.0 1.0 2.0 mesh mesh_ycme 0.0 2.0 4.0
mesh mesh_xfm  10    mesh mesh_yfm  20 # all arrays are read in 1d;
logically these ordered by x->y->z mesh mesh_mat  0 0 1 0     #
Material    material number_materials 2 # moderator    material
sigma_t    0   0.1890 1.4633    material sigma_s    0 0 0.1507 0.0000
material sigma_s    0 1 0.0380 1.4536 # fuel    material sigma_t    1
0.2263 1.0119    material sigma_s 1 0 0.2006 0.0000    material
sigma_s    1 1 0.0161 0.9355    material nu_sigma_f 1   0.0067 0.1241
material chi        1   1.0000 0.0000 # Fixed source (only isotropic
or constant)    source type isotropic source number_spectra 2
source spectrum 0 1.0 0.0    source spectrum 1 0.0 0.0    source map
0 0 1 0   *

See:  InputDB

C++ includes: StupidParser.hh ";

%feature("docstring")  detran::StupidParser::StupidParser "detran::StupidParser::StupidParser(int argc, char **argv)

Constructor.

Parameters:
-----------

argc:  Number of arguments

argv:  Arguments ";

%feature("docstring")  detran::StupidParser::parse_input "StupidParser::SP_input detran::StupidParser::parse_input()

Parse the input. ";

%feature("docstring")  detran::StupidParser::parse_mesh "StupidParser::SP_mesh detran::StupidParser::parse_mesh()

Parse mesh. ";

%feature("docstring")  detran::StupidParser::parse_material "StupidParser::SP_material detran::StupidParser::parse_material()

Parse material. ";


// File: classdetran_1_1Sweeper.xml
%feature("docstring") detran::Sweeper "

Sweeper for discrete ordinates problems.

The within-group transport equation is \\\\[ \\\\mathbf{L}\\\\psi = Q
\\\\, , \\\\] where $ \\\\mathbf{L} $ is the streaming and collision
operator and $ Q $ is a discrete representation of all source
contributions.

To invert the operator $ \\\\mathbf{L} $, we \"sweep\" over the mesh
for all angles, which gives us updated angular fluxes in each cell.
Actually, the flux *moments* are updated, while the discrete angular
flux is optionally stored.

Relevant input database entries: store_angular_flux [int]

equation [string]

C++ includes: Sweeper.hh ";

%feature("docstring")  detran::Sweeper::Sweeper "home robertsj
Research detran source src transport Sweeper cc detran::Sweeper< D
>::Sweeper(SP_input input, SP_mesh mesh, SP_material material,
SP_quadrature quadrature, SP_state state, SP_boundary boundary,
SP_sweepsource sweepsource)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

material:  Material database.

quadrature:  Angular quadrature.

state:   State vectors. ";

%feature("docstring")  detran::Sweeper::~Sweeper "virtual
detran::Sweeper< D >::~Sweeper()

Virtual destructor. ";

%feature("docstring")  detran::Sweeper::sweep "virtual void
detran::Sweeper< D >::sweep(moments_type &phi)
__attribute__((always_inline))=0

Sweep over all angles and space.

Note, if the angular flux is to be updated, it is done directly to via
State. Having sweep take the flux as an explicit argument allows
various input types (e.g. Krylov vectors) without having to go through
State. The angular flux will never be a direct unknown.

Parameters:
-----------

phi:  reference to moments vector to be updated ";

%feature("docstring")  detran::Sweeper::setup_group "void
detran::Sweeper< D >::setup_group(const size_t g)

Setup the equations for the group. Default sets only the group index.
";

%feature("docstring")  detran::Sweeper::set_update_psi "void
detran::Sweeper< D >::set_update_psi(const bool v)

Allows the psi update to occur whenever needed. ";

%feature("docstring")  detran::Sweeper::set_update_boundary "void
detran::Sweeper< D >::set_update_boundary(const bool v)

Switch on-the-fly boundary updates on or off. ";

%feature("docstring")  detran::Sweeper::update_boundary "bool
detran::Sweeper< D >::update_boundary() const ";

%feature("docstring")  detran::Sweeper::number_sweeps "detran_utilities::size_t detran::Sweeper< D >::number_sweeps() const
";

%feature("docstring")  detran::Sweeper::set_adjoint "void
detran::Sweeper< D >::set_adjoint(const bool adjoint)

Set adjoint. ";

%feature("docstring")  detran::Sweeper::is_adjoint "bool
detran::Sweeper< D >::is_adjoint() const

Is adjoint? ";

%feature("docstring")  detran::Sweeper::set_tally "void
detran::Sweeper< D >::set_tally(SP_tally tally)

Set a boundary flux tally. ";


// File: classdetran_1_1Sweeper1D.xml
%feature("docstring") detran::Sweeper1D "

Sweeper for 1D discrete ordinates problems.

C++ includes: Sweeper1D.hh ";

%feature("docstring")  detran::Sweeper1D::Sweeper1D "home robertsj
Research detran source src transport Sweeper1D cc detran::Sweeper1D<
EQ >::Sweeper1D(SP_input input, SP_mesh mesh, SP_material material,
SP_quadrature quadrature, SP_state state, SP_boundary boundary,
SP_sweepsource sweepsource)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

material:  Material database.

quadrature:  Angular quadrature.

state:   State vectors.

boundary:  Boundary based on mesh.

sweepsource:  Sweep source constructor. ";

%feature("docstring")  detran::Sweeper1D::~Sweeper1D "virtual
detran::Sweeper1D< EQ >::~Sweeper1D()

Virtual destructor. ";

%feature("docstring")  detran::Sweeper1D::sweep "void
detran::Sweeper1D< EQ >::sweep(State::moments_type &phi)

Sweep. ";


// File: classdetran_1_1Sweeper2D.xml
%feature("docstring") detran::Sweeper2D "

Sweeper for 2D discrete ordinates problems.

C++ includes: Sweeper2D.hh ";

%feature("docstring")  detran::Sweeper2D::Sweeper2D "home robertsj
Research detran source src transport Sweeper2D cc detran::Sweeper2D<
EQ >::Sweeper2D(SP_input input, SP_mesh mesh, SP_material material,
SP_quadrature quadrature, SP_state state, SP_boundary boundary,
SP_sweepsource sweepsource)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

material:  Material database.

quadrature:  Angular quadrature.

state:   State vectors.

boundary:  Boundary based on mesh.

sweepsource:  Sweep source constructor. ";

%feature("docstring")  detran::Sweeper2D::~Sweeper2D "virtual
detran::Sweeper2D< EQ >::~Sweeper2D()

Virtual destructor. ";

%feature("docstring")  detran::Sweeper2D::sweep "void
detran::Sweeper2D< EQ >::sweep(moments_type &phi)
__attribute__((always_inline))

Sweep. ";


// File: classdetran_1_1Sweeper2DMOC.xml
%feature("docstring") detran::Sweeper2DMOC "

Sweeper for 2D MOC problems.

C++ includes: Sweeper2DMOC.hh ";

%feature("docstring")  detran::Sweeper2DMOC::Sweeper2DMOC "home
robertsj Research detran source src transport Sweeper2DMOC cc
detran::Sweeper2DMOC< EQ >::Sweeper2DMOC(SP_input input, SP_mesh mesh,
SP_material material, SP_quadrature quadrature, SP_state state,
SP_boundary boundary, SP_sweepsource sweepsource)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Tracked mesh.

material:  Material database.

quadrature:  Angular quadrature for MOC.

state:   State vectors.

boundary:  Boundary based on tracks.

sweepsource:  Sweep source constructor. ";

%feature("docstring")  detran::Sweeper2DMOC::~Sweeper2DMOC "virtual
detran::Sweeper2DMOC< EQ >::~Sweeper2DMOC()

Virtual destructor. ";

%feature("docstring")  detran::Sweeper2DMOC::sweep "void
detran::Sweeper2DMOC< EQ >::sweep(moments_type &phi)

Sweep. ";


// File: classdetran_1_1Sweeper3D.xml
%feature("docstring") detran::Sweeper3D "

Sweeper for 3D discrete ordinates problems.

C++ includes: Sweeper3D.hh ";

%feature("docstring")  detran::Sweeper3D::Sweeper3D "home robertsj
Research detran source src transport Sweeper3D cc detran::Sweeper3D<
EQ >::Sweeper3D(SP_input input, SP_mesh mesh, SP_material material,
SP_quadrature quadrature, SP_state state, SP_boundary boundary,
SP_sweepsource sweepsource)

Constructor.

Parameters:
-----------

input:  User input database.

mesh:  Cartesian mesh.

material:  Material database.

quadrature:  Angular quadrature.

state:   State vectors.

boundary:  Boundary based on mesh.

sweepsource:  Sweep source constructor. ";

%feature("docstring")  detran::Sweeper3D::~Sweeper3D "virtual
detran::Sweeper3D< EQ >::~Sweeper3D()

Virtual destructor. ";

%feature("docstring")  detran::Sweeper3D::sweep "void
detran::Sweeper3D< EQ >::sweep(moments_type &phi)

Sweep. ";


// File: classdetran_1_1SweepOperator.xml
%feature("docstring") detran::SweepOperator "

Operator wrapped around a transport sweep.

C++ includes: SweepOperator.hh ";

%feature("docstring")  detran::SweepOperator::SweepOperator "detran::SweepOperator< D >::SweepOperator(SP_state state, SP_boundary
boundary, SP_sweeper sweeper, SP_sweepsource source) ";

%feature("docstring")  detran::SweepOperator::~SweepOperator "virtual
detran::SweepOperator< D >::~SweepOperator() ";

%feature("docstring")  detran::SweepOperator::display "void
detran::SweepOperator< D >::display() const ";

%feature("docstring")  detran::SweepOperator::multiply "void
detran::SweepOperator< D >::multiply(const Vector &x, Vector &y) ";

%feature("docstring")  detran::SweepOperator::multiply_transpose "virtual void detran::SweepOperator< D >::multiply_transpose(const
Vector &x, Vector &y) ";


// File: classdetran_1_1SweepSource.xml
%feature("docstring") detran::SweepSource "

Construct the source for sweeping.

This class defines a general right-hand-side construct for sweeps
(which are the basis for all Detran transport solvers). That is, the
SweepSource is a source for a sweep along a particular angle.

The source for a given sweep is comprised in general of these four
components: within-group scattering source

in-scattering source (from both up- and down-scattering)

fission source

external source Algorithmically, the fission source behaves as an
external source within an inner iteration. Moreover, and of importance
for response function generation, there can be a fission source (i.e.
multiplication) together with an fixed source.

Recall that a within group equation is represented as \\\\[
\\\\mathbf{T}[\\\\psi]_g =
[\\\\mathbf{M}][\\\\mathbf{S}]_{gg}[\\\\phi]_g + \\\\bar{Q}_g \\\\, ,
\\\\] where $ \\\\bar{Q}_g $ represents everything but the within-
group scattering. The right hand side is the sweep source (but not the
right hand side for the linear system of interest, which casts the
problem in terms of flux moments!).

As a sanity check, consider a purely isotropic flux in one group with
isotropic scattering. Then we must have the equation \\\\[
\\\\mathbf{T}\\\\psi = \\\\Sigma_{0} \\\\phi_{00} / 4\\\\pi \\\\, ,
\\\\] which indicates that application of $ \\\\mathbf{M} $ implies a
normalization.

For each inner iteration, the scattering components are stored in
moments form. The M operator is then applied using one row at a time
for the associated angle being swept.

Keep in mind where we're headed: solving the linear system $
(\\\\mathbf{I}-\\\\mathbf{D}\\\\mathbf{T}^{-1}\\\\mathbf{M}\\\\mathb
f{S})\\\\phi = \\\\mathbf{D}\\\\mathbf{T}^{-1}q $. Hence, the sweep
source is q, and we'll invert the transport operator T via sweeps.

See:   ScatterSource, FissionSource, ExternalSource

C++ includes: SweepSource.hh ";

%feature("docstring")  detran::SweepSource::SweepSource "detran::SweepSource< D >::SweepSource(SP_state state, SP_mesh mesh,
SP_quadrature quadrature, SP_material material, SP_MtoD MtoD, bool
implicit_fission=false)

Constructor.

Parameters:
-----------

state:  The state.

mesh:  The mesh.

angularmesh:  The angular mesh.

materials:  The material library.

momentsindex:  The moments index.

m_operator:  The moments to discrete operator. ";

%feature("docstring")  detran::SweepSource::set_moment_source "void
detran::SweepSource< D >::set_moment_source(SP_externalsource source)

Set an external moment source.

Parameters:
-----------

source:  Smart pointer to external source ";

%feature("docstring")  detran::SweepSource::set_discrete_source "void
detran::SweepSource< D >::set_discrete_source(SP_externalsource
source)

Set an external discrete source.

There really is no difference between a moment and discrete source in
representation (the functions are the same), but the sweep source adds
the discrete portion of the source before giving the client the
source. This way, the moment sources can be done once and the discrete
component can be built as needed.

Parameters:
-----------

source:  Smart pointer to external source ";

%feature("docstring")  detran::SweepSource::set_fission_source "void
detran::SweepSource< D >::set_fission_source(SP_fissionsource source)

Set a fission source.

Parameters:
-----------

source:  Smart pointer to fission source ";

%feature("docstring")  detran::SweepSource::get_scatter_source "SP_scattersource detran::SweepSource< D >::get_scatter_source() ";

%feature("docstring")  detran::SweepSource::build_fixed "void
detran::SweepSource< D >::build_fixed(const size_t g)

Add all fixed moments source.

Creates a fixed group source using external moments sources and, if
applicable, fission sources. Note, fission is only fixed within an
eigenproblem. ";

%feature("docstring")  detran::SweepSource::build_fixed_with_scatter "void detran::SweepSource< D >::build_fixed_with_scatter(const size_t
g)

Build fixed with in-scatter.

Creates a fixed group source as build_fixed, but adds the in-scatter
source. This is the typical construction to be called within source
iteration. ";

%feature("docstring")
detran::SweepSource::build_fixed_with_downscatter "void
detran::SweepSource< D >::build_fixed_with_downscatter(const size_t g,
const size_t g_cutoff)

Build fixed with downs-scatter.

Creates a fixed group source as build_fixed, but adds the in-scatter
source. This is the typical construction to be called within source
iteration.

Parameters:
-----------

g:  Group of current solve

g_cutoff:  Last group to include in downscatter source ";

%feature("docstring")  detran::SweepSource::build_within_group_scatter
"void detran::SweepSource< D >::build_within_group_scatter(const
size_t g, const moments_type &phi)

Build within-group scattering source.

Called before each sweep. ";

%feature("docstring")  detran::SweepSource::build_total_scatter "void
detran::SweepSource< D >::build_total_scatter(const size_t g, const
size_t g_cutoff, const State::vec_moments_type &phi)

Build total scattering source.

This builds the complete group scatter source using the given
multigroup flux vector. This routine would find use in a multigroup
Krylov solver. ";

%feature("docstring")  detran::SweepSource::reset "void
detran::SweepSource< D >::reset()

Reset all the internal source vectors to zero. ";

%feature("docstring")  detran::SweepSource::source "void
detran::SweepSource< D >::source(const size_t g, const size_t o, const
size_t a, sweep_source_type &s)

Fill a source vector. ";

%feature("docstring")  detran::SweepSource::fixed_group_source "const
moments_type& detran::SweepSource< D >::fixed_group_source() const

Return the fixed source for the current group. ";

%feature("docstring")  detran::SweepSource::scatter_group_source "const moments_type& detran::SweepSource< D >::scatter_group_source()
const

Return the scatter source for the current group. ";


// File: classdetran_1_1SyntheticDiscreteSource.xml
%feature("docstring") detran::SyntheticDiscreteSource "C++ includes:
SyntheticDiscreteSource.hh ";

%feature("docstring")
detran::SyntheticDiscreteSource::SyntheticDiscreteSource "home
robertsj Research detran source src kinetics SyntheticDiscreteSource
cc detran::SyntheticDiscreteSource::SyntheticDiscreteSource(const
size_t number_groups, SP_mesh mesh, SP_quadrature quadrature,
SP_material material)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Mesh

quadrature:  Quadrature

material:  Material ";

%feature("docstring")  detran::SyntheticDiscreteSource::source "double detran::SyntheticDiscreteSource::source(const size_t cell,
const size_t group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran::SyntheticDiscreteSource::source "double detran::SyntheticDiscreteSource::source(const size_t cell,
const size_t group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";

%feature("docstring")  detran::SyntheticDiscreteSource::build "void
detran::SyntheticDiscreteSource::build(const double dt, const
vec_states &states, const vec_precursors &precursors, const size_t
order)

Build the synthetic source given previous iterates.

Depending on the order of the method, we may use ";


// File: classdetran_1_1SyntheticMomentSource.xml
%feature("docstring") detran::SyntheticMomentSource "C++ includes:
SyntheticMomentSource.hh ";

%feature("docstring")
detran::SyntheticMomentSource::SyntheticMomentSource "home robertsj
Research detran source src kinetics SyntheticMomentSource cc
detran::SyntheticMomentSource::SyntheticMomentSource(const size_t
number_groups, SP_mesh mesh, SP_material material)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Mesh

material:  Material ";

%feature("docstring")  detran::SyntheticMomentSource::source "double
detran::SyntheticMomentSource::source(const size_t cell, const size_t
group)

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran::SyntheticMomentSource::source "double
detran::SyntheticMomentSource::source(const size_t cell, const size_t
group, const size_t angle)

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";

%feature("docstring")  detran::SyntheticMomentSource::build "void
detran::SyntheticMomentSource::build(const double dt, const vec_states
&states, const vec_precursors &precursors, const size_t order=1)

Build the synthetic source given previous iterates.

Depending on the order of the method, we may use ";


// File: classdetran_1_1SyntheticSource.xml
%feature("docstring") detran::SyntheticSource "

Time step source contribution from previous step solution.

In all the time stepping schemes implemented in Detran, the time step
is cast as a fixed source problem with a synthetic source. This class
builds the source from previous state vectors (using either the
angular or scalar flux) and previous precursor vectors.

This class implements synthetic sources based on BDF discretizations
of orders 1-6. Note, first order BDF is equivalent to backward Euler
and is used (via a half-step) to implement the implicit midpoint rule
as used in PARTISN.

C++ includes: SyntheticSource.hh ";

%feature("docstring")  detran::SyntheticSource::SyntheticSource "detran::SyntheticSource::SyntheticSource(const size_t number_groups,
SP_mesh mesh, SP_quadrature quadrature, SP_material material, bool
discrete=false)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Mesh

quadrature:  Angular mesh

material:  Material

discrete:  Flag for discrete source ";

%feature("docstring")  detran::SyntheticSource::build "virtual void
detran::SyntheticSource::build(const double dt, const vec_states
&states, const vec_precursors &precursors, const size_t order=1)=0

Build the synthetic source given previous iterates.

Depending on the order of the method, we may use ";


// File: classdetran__angle_1_1TabuchiYamamoto.xml
%feature("docstring") detran_angle::TabuchiYamamoto "

Quadrature of Tabuchi and Yamamoto for polar integration.

Tabuchi and Yamamoto define an optimal polar quadrature such that the
maximum error of an approximate Ki(3,x) over all x is minimized.
Ki(3,x) is a Bickley-Naylor function, which is related to modified
Bessel functions.

Their original paper gives quadratures for 1, 2, and 3 abscissa. We
extend that to 5 abscissa. Additionally, we implement another size
abscissa based on an L2 minimization of the error (rather than the
infinity norm).

C++ includes: TabuchiYamamoto.hh ";

%feature("docstring")  detran_angle::TabuchiYamamoto::TabuchiYamamoto
"detran_angle::TabuchiYamamoto::TabuchiYamamoto(size_t number_polar)
";

%feature("docstring")  detran_angle::TabuchiYamamoto::~TabuchiYamamoto
"detran_angle::TabuchiYamamoto::~TabuchiYamamoto() ";


// File: classdetran__test_1_1TestDriver.xml
%feature("docstring") detran_test::TestDriver "

Drives a test set.

C++ includes: TestDriver.hh ";


// File: classcallow_1_1TestMatrixShell.xml
%feature("docstring") callow::TestMatrixShell "C++ includes:
matrixshell_fixture.hh ";

%feature("docstring")  callow::TestMatrixShell::TestMatrixShell "callow::TestMatrixShell::TestMatrixShell(const size_t m) ";

%feature("docstring")  callow::TestMatrixShell::multiply "void
callow::TestMatrixShell::multiply(const Vector &x, Vector &y) ";

%feature("docstring")  callow::TestMatrixShell::multiply_transpose "void callow::TestMatrixShell::multiply_transpose(const Vector &x,
Vector &y) ";


// File: classdetran_1_1TimeDependentExternalSource.xml
%feature("docstring") detran::TimeDependentExternalSource "

Base class for time-dependent external sources.

The only addition beyond the base external source class is to specify
a time at which the source will be evaluated during the transport
solve for the current step.

C++ includes: TimeDependentExternalSource.hh ";

%feature("docstring")
detran::TimeDependentExternalSource::TimeDependentExternalSource "det
ran::TimeDependentExternalSource::TimeDependentExternalSource(size_t
number_groups, SP_mesh mesh, SP_quadrature quadrature, bool discrete)

Constructor.

Parameters:
-----------

number_groups:  Number of energy groups

mesh:  Pointer to mesh

quadrature:  Pointer to angular quadrature

discrete:  Flag to indicate treatment as moment or discrete ";

%feature("docstring")
detran::TimeDependentExternalSource::~TimeDependentExternalSource "virtual
detran::TimeDependentExternalSource::~TimeDependentExternalSource()

Virtual destructor. ";

%feature("docstring")  detran::TimeDependentExternalSource::source "virtual double detran::TimeDependentExternalSource::source(const
size_t cell, const size_t group)=0

Get moments source for cell.

Units are n/cc-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group ";

%feature("docstring")  detran::TimeDependentExternalSource::source "virtual double detran::TimeDependentExternalSource::source(const
size_t cell, const size_t group, const size_t angle)=0

Get discrete source for cell and cardinal angle.

Units are n/cc-ster-sec

Parameters:
-----------

cell:  Mesh cell

group:  Energy group

angle:  Cardinal angle index ";

%feature("docstring")  detran::TimeDependentExternalSource::set_time "virtual void detran::TimeDependentExternalSource::set_time(const
double time)=0

Set the time of the source.

Parameters:
-----------

time:  Time at which source is evaluated ";


// File: classdetran_1_1TimeDependentManager.xml
%feature("docstring") detran::TimeDependentManager "

Manage solution of time-dependnent problems.

C++ includes: TimeDependentManager.hh ";

%feature("docstring")
detran::TimeDependentManager::TimeDependentManager "home robertsj
Research detran source src solvers TimeDependentManager cc
detran::TimeDependentManager< D >::TimeDependentManager(SP_input
input, SP_material material, SP_mesh mesh)

Constructor.

Parameters:
-----------

input:  parameter database

material:  material database

mesh:  mesh definition ";

%feature("docstring")  detran::TimeDependentManager::setup "void
detran::TimeDependentManager< D >::setup()

Sets up a problem to be solved.

This sets up the discretization, quadrature, boundary, and state. By
changing the appropriate database parameters, a call to setup will
rebuild the problem using a different discretization, etc.

Source defining external sources, the quadrature is often needed. By
calling setup first, the client can extract the quadrature for
building an external source. ";

%feature("docstring")  detran::TimeDependentManager::set_source "void
detran::TimeDependentManager< D >::set_source(SP_source q)

Add a new external source.

This is called after setup, since external sources often need the
quadrature set produced by the manager. ";

%feature("docstring")  detran::TimeDependentManager::set_solver "bool
detran::TimeDependentManager< D >::set_solver()

Set the solver based on parameter database.

This must be called after setup is called and after all sources have
been set. By changing the appropriate parameters, a problem can be re-
run with a different solver. ";

%feature("docstring")  detran::TimeDependentManager::solve "bool
detran::TimeDependentManager< D >::solve(const double keff=1.0)

Solve the system.

Parameters:
-----------

keff:  Scaling factor for multiplying problems ";


// File: classdetran_1_1TimeDependentMaterial.xml
%feature("docstring") detran::TimeDependentMaterial "

Base class for time-dependent materials.

This implementation allows for arbitrary time-dependent materials to
be defined, possibly as a function of state.

The client must implement the private implementation of update, which
defines the normal cross section values. This private function is
called in update, which then adjusts the cross sections for the time
step (i.e. it produces synthetic cross sections that are dependent on
the time increment)

C++ includes: TimeDependentMaterial.hh ";

%feature("docstring")
detran::TimeDependentMaterial::TimeDependentMaterial "detran::TimeDependentMaterial::TimeDependentMaterial(const size_t
number_materials, const size_t number_energy_groups, const size_t
number_precursor_groups, std::string name=\"TimeDependentMaterial\")

Constructor.

Parameters:
-----------

number_materials:  Number of materials

number_energy_groups:  Number of energy groups

number_precursor_groups:  Number of precursor groups

state:   State vector ";

%feature("docstring")
detran::TimeDependentMaterial::~TimeDependentMaterial "virtual
detran::TimeDependentMaterial::~TimeDependentMaterial()

Virtual destructor. ";

%feature("docstring")  detran::TimeDependentMaterial::set_eigenvalue "void detran::TimeDependentMaterial::set_eigenvalue(const double keff)
";

%feature("docstring")  detran::TimeDependentMaterial::set_state "void
detran::TimeDependentMaterial::set_state(SP_state s)

Set the state vector. ";

%feature("docstring")  detran::TimeDependentMaterial::state "SP_state
detran::TimeDependentMaterial::state()

Get the state vector. ";

%feature("docstring")  detran::TimeDependentMaterial::time "double
detran::TimeDependentMaterial::time() const ";

%feature("docstring")  detran::TimeDependentMaterial::dt "double
detran::TimeDependentMaterial::dt() const ";

%feature("docstring")  detran::TimeDependentMaterial::order "size_t
detran::TimeDependentMaterial::order() const ";

%feature("docstring")  detran::TimeDependentMaterial::kcrit "double
detran::TimeDependentMaterial::kcrit() const ";

%feature("docstring")  detran::TimeDependentMaterial::update "void
detran::TimeDependentMaterial::update(const double t, const double dt,
const size_t order, const bool flag=true)

Update the materials.

Parameters:
-----------

t:  New time (in seconds)

dt:  Time step (in seconds)

order:  BDF order  synthetic Flag for creating a synthetic material ";


// File: classdetran__utilities_1_1Timer.xml
%feature("docstring") detran_utilities::Timer "C++ includes: Timer.hh
";

%feature("docstring")  detran_utilities::Timer::Timer "detran_utilities::Timer::Timer() ";

%feature("docstring")  detran_utilities::Timer::tic "void
detran_utilities::Timer::tic()

Begin the timer around a single code block. ";

%feature("docstring")  detran_utilities::Timer::toc "double
detran_utilities::Timer::toc(bool flag=false)

Return the time elapsed around a single code block. ";

%feature("docstring")  detran_utilities::Timer::wtime "double
detran_utilities::Timer::wtime()

Return the current wall time (only differences are meaningful) ";

%feature("docstring")  detran_utilities::Timer::function_tic "void
detran_utilities::Timer::function_tic()

Begin the timer at the beginning of a function call. ";

%feature("docstring")  detran_utilities::Timer::function_toc "void
detran_utilities::Timer::function_toc(std::string function)

Log the function time. ";


// File: classdetran_1_1TimeStepper.xml
%feature("docstring") detran::TimeStepper "

Solve a time dependent problem by stepping forward in time.

Available time stepping options are the BDF methods of order 1 through
6 and the implicit midpoint rule.

C++ includes: TimeStepper.hh ";

%feature("docstring")  detran::TimeStepper::TimeStepper "detran::TimeStepper< D >::TimeStepper(SP_input input, SP_material
material, SP_mesh mesh, bool multiply)

Constructor.

Parameters:
-----------

input:

material:  ";

%feature("docstring")  detran::TimeStepper::~TimeStepper "virtual
detran::TimeStepper< D >::~TimeStepper()

Virtual destructor. ";

%feature("docstring")  detran::TimeStepper::add_source "void
detran::TimeStepper< D >::add_source(SP_tdsource source)

Add a time-dependent external source. ";

%feature("docstring")  detran::TimeStepper::solve "void
detran::TimeStepper< D >::solve(SP_state initial_state)

Solve. ";

%feature("docstring")  detran::TimeStepper::state "SP_state
detran::TimeStepper< D >::state()

Getters. ";

%feature("docstring")  detran::TimeStepper::mesh "SP_mesh
detran::TimeStepper< D >::mesh() ";

%feature("docstring")  detran::TimeStepper::material "SP_material
detran::TimeStepper< D >::material() ";

%feature("docstring")  detran::TimeStepper::quadrature "SP_quadrature
detran::TimeStepper< D >::quadrature() ";

%feature("docstring")  detran::TimeStepper::monitor_level "size_t
detran::TimeStepper< D >::monitor_level() const ";

%feature("docstring")  detran::TimeStepper::precursor "SP_precursors
detran::TimeStepper< D >::precursor() ";

%feature("docstring")  detran::TimeStepper::multiphysics "SP_multiphysics detran::TimeStepper< D >::multiphysics() ";

%feature("docstring")  detran::TimeStepper::fissionsource "SP_fissionsource detran::TimeStepper< D >::fissionsource() ";

%feature("docstring")  detran::TimeStepper::residual_norm "double
detran::TimeStepper< D >::residual_norm() ";

%feature("docstring")  detran::TimeStepper::set_monitor "void
detran::TimeStepper< D >::set_monitor(monitor_pointer monitor, void
*monitor_data=NULL)

Set a user-defined monitor function. ";

%feature("docstring")  detran::TimeStepper::set_multiphysics "void
detran::TimeStepper< D >::set_multiphysics(SP_multiphysics ic,
multiphysics_pointer update_multiphysics_rhs, void
*multiphysics_data=NULL)

Set the multiphysics.

The user must initialize the multiphysics state vector. The user also
specifies the multiphysics update function as well as any data
required. ";


// File: classdetran__geometry_1_1Track.xml
%feature("docstring") detran_geometry::Track "

Represents a track across a domain, consisting of several segments.

C++ includes: Track.hh ";

%feature("docstring")  detran_geometry::Track::Track "detran_geometry::Track::Track(Point r0, Point r1)

Constructor.

Parameters:
-----------

r0:  Entrance point

r1:  Exit point ";

%feature("docstring")  detran_geometry::Track::enter "Point
detran_geometry::Track::enter() const

Return my entrance point. ";

%feature("docstring")  detran_geometry::Track::exit "Point
detran_geometry::Track::exit() const

Return my exit point. ";

%feature("docstring")  detran_geometry::Track::number_segments "int
detran_geometry::Track::number_segments() const ";

%feature("docstring")  detran_geometry::Track::add_segment "void
detran_geometry::Track::add_segment(Segment s) ";

%feature("docstring")  detran_geometry::Track::segment "Segment&
detran_geometry::Track::segment(size_t i)

Mutable reference to segment. ";

%feature("docstring")  detran_geometry::Track::segment "const
Segment& detran_geometry::Track::segment(size_t i) const

Const reference to segment. ";

%feature("docstring")  detran_geometry::Track::begin "iterator
detran_geometry::Track::begin()

Iterator to the beginning of the track. ";

%feature("docstring")  detran_geometry::Track::rbegin "riterator
detran_geometry::Track::rbegin()

Iterator to the end of the track. ";

%feature("docstring")  detran_geometry::Track::cos_phi "double
detran_geometry::Track::cos_phi() const

Return the track cosine (with respect to x) ";

%feature("docstring")  detran_geometry::Track::sin_phi "double
detran_geometry::Track::sin_phi() const

Return the track sine (with respect to x) ";


// File: classdetran__geometry_1_1TrackDB.xml
%feature("docstring") detran_geometry::TrackDB "

Database of tracks.

The track database contains all the information needed to represent
the MOC problem with respect to space and angle.

Currently, this is limited to 2-D tracking. Tracks are only assigned
for azimuths over [0, pi], using symmetry for [pi, 2*pi].

Tracks are stored in order in order of decreasing y and increasing x
for 0, pi/2 and decreasing x for pi/2, pi. This is how tracks are
indexed. The track index in the other two octancts keeps the index of
their reflection.

C++ includes: TrackDB.hh ";

%feature("docstring")  detran_geometry::TrackDB::TrackDB "detran_geometry::TrackDB::TrackDB(int num_azimuths, int num_regions,
SP_quadrature quad)

Constructor.

Parameters:
-----------

num_azimuths:  Number of azimuths in first two octants ";

%feature("docstring")  detran_geometry::TrackDB::~TrackDB "detran_geometry::TrackDB::~TrackDB() ";

%feature("docstring")  detran_geometry::TrackDB::add_track "void
detran_geometry::TrackDB::add_track(size_t a, SP_track t) ";

%feature("docstring")  detran_geometry::TrackDB::setup_angle "void
detran_geometry::TrackDB::setup_angle(size_t a, double c_phi, double
s_phi, double space) ";

%feature("docstring")  detran_geometry::TrackDB::normalize "void
detran_geometry::TrackDB::normalize(vec_dbl &volume)

Normalize the tracks given a vector of true volumes. ";

%feature("docstring")  detran_geometry::TrackDB::display "void
detran_geometry::TrackDB::display() const

Pretty display of all track. ";


// File: classdetran__geometry_1_1Tracker.xml
%feature("docstring") detran_geometry::Tracker "

Track a mesh.

C++ includes: Tracker.hh ";

%feature("docstring")  detran_geometry::Tracker::Tracker "detran_geometry::Tracker::Tracker(SP_mesh mesh, SP_quadrature
quadrature) ";

%feature("docstring")  detran_geometry::Tracker::trackdb "SP_trackdb
detran_geometry::Tracker::trackdb() const ";

%feature("docstring")  detran_geometry::Tracker::meshmoc "SP_mesh
detran_geometry::Tracker::meshmoc() ";

%feature("docstring")  detran_geometry::Tracker::normalize "void
detran_geometry::Tracker::normalize() ";


// File: classdetran_1_1TransportManager.xml
%feature("docstring") detran::TransportManager "

Base class for all transport managers.

The purpose of an untemplated base mesh is to simplify initialization
of libraries, etc.

C++ includes: TransportManager.hh ";

%feature("docstring")  detran::TransportManager::TransportManager "detran::TransportManager::TransportManager(int argc, char *argv[])

Constructor.

Initializes external libraries

Parameters:
-----------

argc:  command line count

argv:  command line values ";

%feature("docstring")  detran::TransportManager::TransportManager "detran::TransportManager::TransportManager()

Default constructor. ";

%feature("docstring")  detran::TransportManager::~TransportManager "virtual detran::TransportManager::~TransportManager()

Destructor. This initializes all external libraries. ";


// File: structcallow_1_1triplet.xml
%feature("docstring") callow::triplet "

Structure for storing COO matrix.

C++ includes: Triplet.hh ";

%feature("docstring")  callow::triplet::triplet "callow::triplet::triplet(int ii=-1, int jj=-1, double vv=0.0) ";


// File: classTWIGLMaterial.xml
%feature("docstring") TWIGLMaterial "";

%feature("docstring")  TWIGLMaterial::TWIGLMaterial "TWIGLMaterial::TWIGLMaterial(int perturbation=0, bool transport=false)
";

%feature("docstring")  TWIGLMaterial::update_impl "void
TWIGLMaterial::update_impl() ";


// File: classdetran__angle_1_1Uniform.xml
%feature("docstring") detran_angle::Uniform "

Uniformly-spaced azimuthal quadrature set.

C++ includes: Uniform.hh ";

%feature("docstring")  detran_angle::Uniform::Uniform "detran_angle::Uniform::Uniform(size_t dim, size_t num_azimuths_octant,
size_t num_space, size_t num_polar, std::string polar)

Constructor.

Parameters:
-----------

dim:  Problem dimension (only 2 supported for now)

num_azimuths_octant:  Number of azimuths per octant

num_space:  Number of tracks per azimuth

num_polar:  Number of polar angles in half space

polar:  Polar quadrature string identifier ";

%feature("docstring")  detran_angle::Uniform::~Uniform "detran_angle::Uniform::~Uniform() ";


// File: classdetran__angle_1_1UniformEqual.xml
%feature("docstring") detran_angle::UniformEqual "

2D/3D Uniform, Equal Weight (UEn) quadrature class.

As mentioned in LevelSymmetric, a fundamental problem inherent to LQn
quadrature is presence of negative weights at high order. These
weights produce unphysical solutions (and may inhibit convergence).
For problems where an increased quadrature order (i.e. more angles) is
required to study convergence or simply to get better answers, we
require an arbitrarily high order, positive weight quadrature.

Here, we implement the \"uniform, equal weight\" (UEn) quadrature of
Carew and Zamonsky. The basic idea is to choose uniform azimuthal
divisions and uniform polar cosines. The result is a product
quadrature that can be extended on-the-fly to arbitrary numbers of
angles.

Here, the quadrature order defines the number of polar angles. We take
the number of azimuthal angles to be twice this number. Hence, the
total number of angles is twice the number of polar angles squared.

Note, the angles are stored with polar as the inner index. This is
requested for all product quadratures.

Carew, J. and Zamonsky. G., Nuclear Science and Engineering 131,
199-207 (1999).

C++ includes: UniformEqual.hh ";

%feature("docstring")  detran_angle::UniformEqual::UniformEqual "detran_angle::UniformEqual::UniformEqual(size_t order, size_t dim)

Constructor.

Parameters:
-----------

order:   Quadrature order.

dim:  Problem dimension ";


// File: classdetran_1_1Vacuum.xml
%feature("docstring") detran::Vacuum "

Vacuum boundary condition for SN calculations.

C++ includes: Vacuum.hh ";

%feature("docstring")  detran::Vacuum::Vacuum "detran::Vacuum< D
>::Vacuum(BoundarySN< D > &boundary, const size_t side, SP_input
input, SP_mesh mesh, SP_quadrature quadrature) ";

%feature("docstring")  detran::Vacuum::set "void detran::Vacuum< D
>::set(const size_t g)

Set initial and/or fixed boundary condition. Vacuum does nothing. ";

%feature("docstring")  detran::Vacuum::update "void detran::Vacuum< D
>::update(const size_t g)

Update a boundary following a sweep. Vacuum does nothing. ";

%feature("docstring")  detran::Vacuum::update "void detran::Vacuum< D
>::update(const size_t g, const size_t o, const size_t a)

Update a boundary for a given angle following a sweep. Vacuum does
nothing. ";


// File: classdetran_1_1VacuumMOC.xml
%feature("docstring") detran::VacuumMOC "

Vacuum boundary condition for MOC.

C++ includes: VacuumMOC.hh ";

%feature("docstring")  detran::VacuumMOC::VacuumMOC "detran::VacuumMOC< D >::VacuumMOC(BoundaryMOC< D > &boundary, const
size_t side, SP_input input, SP_mesh mesh, SP_quadrature quadrature)
";

%feature("docstring")  detran::VacuumMOC::set "void
detran::VacuumMOC< D >::set(const size_t g)

Set initial and/or fixed boundary condition. Vacuum does nothing. ";

%feature("docstring")  detran::VacuumMOC::update "void
detran::VacuumMOC< D >::update(const size_t g)

Update a boundary following a sweep. Vacuum does nothing. ";

%feature("docstring")  detran::VacuumMOC::update "void
detran::VacuumMOC< D >::update(const size_t g, const size_t o, const
size_t a)

Update a boundary for a given angle following a sweep. Vacuum does
nothing. ";


// File: structdetran__utilities_1_1vec1__T.xml
%feature("docstring") detran_utilities::vec1_T "C++ includes:
Definitions.hh ";


// File: structdetran__utilities_1_1vec2__T.xml
%feature("docstring") detran_utilities::vec2_T "C++ includes:
Definitions.hh ";


// File: structdetran__utilities_1_1vec3__T.xml
%feature("docstring") detran_utilities::vec3_T "C++ includes:
Definitions.hh ";


// File: structdetran__utilities_1_1vec4__T.xml
%feature("docstring") detran_utilities::vec4_T "C++ includes:
Definitions.hh ";


// File: classcallow_1_1Vector.xml
%feature("docstring") callow::Vector "

Dense vector object.

C++ includes: Vector.hh ";

%feature("docstring")  callow::Vector::~Vector "callow::Vector::~Vector()

Virtual destructor. ";

%feature("docstring")  callow::Vector::resize "void
callow::Vector::resize(const int n, const double v=0.0)

Wipe out the contents and resize. ";

%feature("docstring")  callow::Vector::value "const double &
callow::Vector::value(const int i) const ";

%feature("docstring")  callow::Vector::value "double &
callow::Vector::value(const int i) ";

%feature("docstring")  callow::Vector::dot "double
callow::Vector::dot(const Vector &x)

Inner product of this vector with vector x. ";

%feature("docstring")  callow::Vector::dot "double
callow::Vector::dot(SP_vector x) ";

%feature("docstring")  callow::Vector::norm "double
callow::Vector::norm(const int type=L2)

Norm of this vector. ";

%feature("docstring")  callow::Vector::norm_residual "double
callow::Vector::norm_residual(const Vector &x, const int type=L2)

Norm of the difference of this and another.

Relative norms are with respect to this vector, and zeros are not
checked. ";

%feature("docstring")  callow::Vector::norm_residual "double
callow::Vector::norm_residual(SP_vector x, const int type=L2) ";

%feature("docstring")  callow::Vector::set "void
callow::Vector::set(const double v)

Set all elements of this vector to a value v. ";

%feature("docstring")  callow::Vector::scale "void
callow::Vector::scale(const double v)

Multiply all elements of this vector by a value v. ";

%feature("docstring")  callow::Vector::add "void
callow::Vector::add(const Vector &x)

Add a vector x to this vector. ";

%feature("docstring")  callow::Vector::add "void
callow::Vector::add(SP_vector x) ";

%feature("docstring")  callow::Vector::subtract "void
callow::Vector::subtract(const Vector &x)

Subtract a vector x from this vector. ";

%feature("docstring")  callow::Vector::subtract "void
callow::Vector::subtract(SP_vector x) ";

%feature("docstring")  callow::Vector::multiply "void
callow::Vector::multiply(const Vector &x)

Multiply this vector pointwise with a vector x. ";

%feature("docstring")  callow::Vector::multiply "void
callow::Vector::multiply(SP_vector x) ";

%feature("docstring")  callow::Vector::divide "void
callow::Vector::divide(const Vector &x)

Multiply this vector pointwise with a vector x. ";

%feature("docstring")  callow::Vector::divide "void
callow::Vector::divide(SP_vector x) ";

%feature("docstring")  callow::Vector::copy "void
callow::Vector::copy(const Vector &x)

Copy a vector x to this vector. ";

%feature("docstring")  callow::Vector::copy "void
callow::Vector::copy(SP_vector x) ";

%feature("docstring")  callow::Vector::add_a_times_x "void
callow::Vector::add_a_times_x(const double a, const Vector &x)

Add a vector x times a scalar a to this vector. ";

%feature("docstring")  callow::Vector::add_a_times_x "void
callow::Vector::add_a_times_x(const double a, SP_vector x) ";

%feature("docstring")  callow::Vector::size "int
callow::Vector::size() const ";

%feature("docstring")  callow::Vector::display "void
callow::Vector::display() const

Pretty print to stdout. ";

%feature("docstring")  callow::Vector::print_matlab "void
callow::Vector::print_matlab(std::string filename=\"vector.out\")
const

Formatted write to ascii for matlab. ";


// File: classdetran_1_1WGDiffusionLossOperator.xml
%feature("docstring") detran::WGDiffusionLossOperator "C++ includes:
WGDiffusionLossOperator.hh ";

%feature("docstring")
detran::WGDiffusionLossOperator::WGDiffusionLossOperator "detran::WGDiffusionLossOperator::WGDiffusionLossOperator(SP_input
input, SP_material material, SP_mesh mesh, size_t group)

Constructor.

Parameters:
-----------

input:  Pointer to input parameters

material:  Pointer to materials

mesh:  Pointer to mesh

group:  Group of operator ";

%feature("docstring")  detran::WGDiffusionLossOperator::construct "void detran::WGDiffusionLossOperator::construct()

Rebuild the matrix based on the present material definitions. ";


// File: classdetran_1_1WGPreconditioner.xml
%feature("docstring") detran::WGPreconditioner "

Preconditioner for within-group equation.

The diffusion operator only operates on the scalar flux, and so
addition of higher order moments will require restriction and
projection operations.

C++ includes: WGPreconditioner.hh ";

%feature("docstring")  detran::WGPreconditioner::WGPreconditioner "detran::WGPreconditioner::WGPreconditioner(SP_input input, SP_material
material, SP_mesh mesh, std::string name =\"WG-PC\")

Constructor.

Parameters:
-----------

input:  Input database

material:  Material database

mesh:  Cartesian mesh ";

%feature("docstring")  detran::WGPreconditioner::~WGPreconditioner "virtual detran::WGPreconditioner::~WGPreconditioner()

virtual destructor ";

%feature("docstring")  detran::WGPreconditioner::set_group "void
detran::WGPreconditioner::set_group(const size_t group)

Set the group for this solve. ";

%feature("docstring")  detran::WGPreconditioner::apply "virtual void
detran::WGPreconditioner::apply(Vector &b, Vector &x)=0

Solve Px = b. ";


// File: classdetran_1_1WGSolver.xml
%feature("docstring") detran::WGSolver "

Solve the within-group transport equation.

The within-group transport equation in operator form is \\\\[
@mathbf{L}\\\\psi = @mathbf{MS}\\\\phi + Q \\\\] where $ \\\\mathbf{L}
$ is the streaming and collision operator, {M} is the moment-to-
discrete operator, {S} is the scattering operator, and Q represents
any source considered fixed, which includes in-scatter, fission, and
external sources.

What we are really after is the scalar flux and possibly its higher
order moments. Consequently, we are able to solve a somewhat different
problem then the within group transport equation above. Let us operate
on both sides by $\\\\mathbf{L}^{-1}$ followed by $ \\\\mathbf{D}$ to
get \\\\[ (\\\\mathbf{I} -
\\\\mathbf{D}\\\\mathbf{L}^{-1}\\\\mathbf{MS})\\\\phi = \\\\mathbf{D}
\\\\mathbf{L}^{-1} Q \\\\, . \\\\] Here, $\\\\mathbf{D}$ is the
discrete-to-moment operator, defined such that $ \\\\phi =
\\\\mathbf{D}\\\\psi $.

Notice this is nothing but a linear system of the form $
\\\\mathbf{A}x = b $ where \\\\[ @mathbf{A} = (\\\\mathbf{I} -
\\\\mathbf{D}\\\\mathbf{L}^{-1}\\\\mathbf{MS}) \\\\] and \\\\[ b =
\\\\mathbf{D} \\\\mathbf{L}^{-1} Q \\\\, . \\\\] Moreover, $ b$ is
just the uncollided flux.

A nice overview of approaches for this inner iteration is given by
Larsen and Morel in  Nuclear Computational Science.

Input parameters specific to WGSolver and derived classes:
inner_max_iters [100]

inner_tolerance [1e-5]

inner_print_out [2], 0=never, 1=final, 2=every interval

inner_print_interval [10]

See:   WGSolverSI, WGSolverGMRES

C++ includes: WGSolver.hh ";

%feature("docstring")  detran::WGSolver::WGSolver "home robertsj
Research detran source src solvers wg WGSolver cc detran::WGSolver< D
>::WGSolver(SP_state state, SP_material material, SP_quadrature
quadrature, SP_boundary boundary, const vec_externalsource &q_e,
SP_fissionsource q_f, bool multiply)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

mat:  Material definitions.

quadrature:  Angular mesh.

boundary:  Boundary fluxes.

external_source:  User-defined external source.

fission_source:  Fission source. ";

%feature("docstring")  detran::WGSolver::~WGSolver "virtual
detran::WGSolver< D >::~WGSolver()

Virtual destructor. ";

%feature("docstring")  detran::WGSolver::solve "virtual void
detran::WGSolver< D >::solve(const size_t g)=0

Solve the within group equation for group g. ";

%feature("docstring")  detran::WGSolver::get_sweepsource "SP_sweepsource detran::WGSolver< D >::get_sweepsource() const ";

%feature("docstring")  detran::WGSolver::get_sweeper "SP_sweeper
detran::WGSolver< D >::get_sweeper() const ";

%feature("docstring")  detran::WGSolver::group "int detran::WGSolver<
D >::group() const ";


// File: classdetran_1_1WGSolverGMRES.xml
%feature("docstring") detran::WGSolverGMRES "

Solve the within-group problem with GMRES.

From WGSolver, we know the within-group problem can be written in
operator notation as \\\\[ (\\\\mathbf{I} -
\\\\mathbf{D}\\\\mathbf{L}^{-1}\\\\mathbf{MS})\\\\phi = \\\\mathbf{D}
\\\\mathbf{L}^{-1} Q \\\\, , \\\\] or \\\\[ \\\\mathbf{A}x = b \\\\, .
\\\\]

This class couples with callow to make available its set of applicable
solvers, the default being GMRES. Other solvers are selected by the
parameter database for callow solvers. PETSc solvers are available
through the callow interface as well.

For Krylov iterations to perform successfully, preconditioning is
often required. A good preconditioner $ {M} $ is in some way
\"similar\" to the operator $\\\\f \\\\mathbf{A} $, and applying its
inverse $ {M}^{-1} $ can be done cheaply.

C++ includes: WGSolverGMRES.hh ";

%feature("docstring")  detran::WGSolverGMRES::WGSolverGMRES "home
robertsj Research detran source src solvers wg WGSolverGMRES cc
detran::WGSolverGMRES< D >::WGSolverGMRES(SP_state state, SP_material
material, SP_quadrature quadrature, SP_boundary boundary, const
vec_externalsource &q_e, SP_fissionsource q_f, bool multiply)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

mat:  Material definitions.

quadrature:  Angular mesh.

boundary:  Boundary fluxes.

external_source:  User-defined external source.

fission_source:  Fission source.

multiply:  Flag for fixed source multiplying problem ";

%feature("docstring")  detran::WGSolverGMRES::solve "void
detran::WGSolverGMRES< D >::solve(const size_t g)

Solve the within group equation. ";


// File: classdetran_1_1WGSolverSI.xml
%feature("docstring") detran::WGSolverSI "

Solve the within-group transport equation via source iteration.

This is a \"hand-coded\" implementation. Essentially the same solver
can be had via callow's Richardson iteration (and via callow's
interface to PETSc's Richardson).

C++ includes: WGSolverSI.hh ";

%feature("docstring")  detran::WGSolverSI::WGSolverSI "home robertsj
Research detran source src solvers wg WGSolverSI cc
detran::WGSolverSI< D >::WGSolverSI(SP_state state, SP_material
material, SP_quadrature quadrature, SP_boundary boundary, const
vec_externalsource &q_e, SP_fissionsource q_f, bool multiply)

Constructor.

Parameters:
-----------

state:   State vectors, etc.

mat:  Material definitions.

quadrature:  Angular mesh.

boundary:  Boundary fluxes.

external_source:  User-defined external source.

fission_source:  Fission source.

multiply:  Flag for fixed source multiplying problem ";

%feature("docstring")  detran::WGSolverSI::solve "void
detran::WGSolverSI< D >::solve(const size_t g)

Solve the within group equation. ";


// File: classdetran_1_1WGTransportOperator.xml
%feature("docstring") detran::WGTransportOperator "

Within-group transport operator.

The within-group transport operator is the the left hand side operator
in the equation for the scalar flux moments.

C++ includes: WGTransportOperator.hh ";

%feature("docstring")
detran::WGTransportOperator::WGTransportOperator "detran::WGTransportOperator< D >::WGTransportOperator(SP_state state,
SP_boundary boundary, SP_sweeper sweeper, SP_sweepsource source) ";

%feature("docstring")
detran::WGTransportOperator::~WGTransportOperator "virtual
detran::WGTransportOperator< D >::~WGTransportOperator() ";

%feature("docstring")  detran::WGTransportOperator::set_group "void
detran::WGTransportOperator< D >::set_group(const size_t g) ";

%feature("docstring")  detran::WGTransportOperator::display "void
detran::WGTransportOperator< D >::display() const ";

%feature("docstring")  detran::WGTransportOperator::multiply "void
detran::WGTransportOperator< D >::multiply(const Vector &x, Vector &y)
";

%feature("docstring")  detran::WGTransportOperator::multiply_transpose
"virtual void detran::WGTransportOperator< D
>::multiply_transpose(const Vector &x, Vector &y) ";


// File: classdetran_1_1WithinGroupAcceleration.xml
%feature("docstring") detran::WithinGroupAcceleration "

Base class for within-group coarse mesh acceleration schemes.

C++ includes: WithinGroupAcceleration.hh ";

%feature("docstring")
detran::WithinGroupAcceleration::WithinGroupAcceleration "detran::WithinGroupAcceleration< D >::WithinGroupAcceleration(SP_input
input, SP_material material, SP_coarsemesh coarsemesh, SP_currenttally
currenttally)

Constructor.

Parameters:
-----------

input:  Input database

material:  Material database

coarsemesh:  Coarse mesh

currenttally:  Current tally ";

%feature("docstring")
detran::WithinGroupAcceleration::~WithinGroupAcceleration "virtual
detran::WithinGroupAcceleration< D >::~WithinGroupAcceleration()

Virtual destructor. ";


// File: namespaceangle.xml


// File: namespaceboundary.xml


// File: namespacecallow.xml
%feature("docstring")  callow::compare_triplet "bool
callow::compare_triplet(const triplet &x, const triplet &y) ";

%feature("docstring")  callow::petsc_ksp_monitor "PetscErrorCode
callow::petsc_ksp_monitor(KSP ksp, PetscInt it, PetscReal rnorm, void
*ctx) ";

%feature("docstring")  callow::test_matrix_1 "Matrix::SP_matrix
callow::test_matrix_1(int n=5) ";

%feature("docstring")  callow::test_matrix_2 "Matrix::SP_matrix
callow::test_matrix_2(int n=10) ";


// File: namespacedetran.xml
%feature("docstring")  detran::apply_inv_P_MG "PetscErrorCode
detran::apply_inv_P_MG(PC pc, Vec x, Vec y)

Apply the PC. ";

%feature("docstring")  detran::ts_default_monitor "void
detran::ts_default_monitor(void *data, TimeStepper< D > *ts, int step,
double t, double dt, int it, bool converged)

Default time step monitor.

Parameters:
-----------

data:  Pointer to arbitrary user data

ts:   TimeStepper pointer

step:  Current time step

t:  Time at step

dt:  Step size

it:  Iteration count (for nonlinear problems) ";


// File: namespacedetran__angle.xml
%feature("docstring")  detran_angle::generate_gc_parameters "void
detran_angle::generate_gc_parameters(detran_utilities::size_t m,
detran_utilities::vec_dbl &x, detran_utilities::vec_dbl &w, bool
normalize=false)

Generate Gauss-Chebyshev parameters.

The Gauss-Chebyshev quadrature approximates \\\\[ \\\\int^{1}_{-1}
\\\\frac{f(x)}{\\\\sqrt{1-x^2}} dx \\\\approx \\\\sum^m_{i=1} W_m
f(x_i) \\\\, . \\\\]

Now, the integral we want is not actually weighted, and so we set the
true weights to be \\\\[ w_i = W_i \\\\sqrt{1-x_i^2} \\\\, . \\\\]

The abscissa of the m-point quadrature are the zeros of the Chebyshev
polynomial of degree m + 1.

However, these weights do not sum to two, as we'd expect when
integrating $ f(x) = 1 $, though the sum does approach two for high
order. Optionally, the user can choose to normalize the weights to
two.

Reference: Hildebrand, Introduction to Numerical Analysis

Parameters:
-----------

m:  number of points (i.e. the quadrature order)

x:  temporary array for abscissa

w:  temporary array for weights ";

%feature("docstring")  detran_angle::generate_gl_parameters "void
detran_angle::generate_gl_parameters(detran_utilities::size_t m,
detran_utilities::vec_dbl &x, detran_utilities::vec_dbl &w)

Generate Gauss-Legendre parameters.

The Gauss-Legendre quadrature approximates \\\\[ \\\\int^{1}_{-1} f(x)
dx \\\\approx \\\\sum^m_{i=1} W_m f(x_i) \\\\, . \\\\]

The abscissa $ x_i $ of the m-point quadrature turn out to be the
zeroes of the Legendre polynomial of degree $ m + 1 $. The weights can
be found numerically, as done in this implementation, which is a
modified version of the function gauleg from  Numerical Recipes in C
by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and & B.P. Flannery,
(Cambridge Univ. Press)

Parameters:
-----------

m:  number of points (i.e. the quadrature order)

x:  temporary array for abscissa

w:  temporary array for weights ";


// File: namespacedetran__external__source.xml


// File: namespacedetran__geometry.xml


// File: namespacedetran__ioutils.xml
%feature("docstring")  detran_ioutils::set_value "void
detran_ioutils::set_value(const T &value_in, typename
compound_type_traits< T >::value_type &value_out)

Set the value to be sent to file given the value in memory.

Parameters:
-----------

value_in:  The value as stored in memory

value_out:  The value to be written to file ";

%feature("docstring")  detran_ioutils::set_value "void
detran_ioutils::set_value(const std::string &value_in,
compound_type_traits< std::string >::value_type &value_out) ";

%feature("docstring")  detran_ioutils::set_value "void
detran_ioutils::set_value(const std::vector< int > &value_in,
compound_type_traits< std::vector< int > >::value_type &value_out) ";

%feature("docstring")  detran_ioutils::set_value "void
detran_ioutils::set_value(const std::vector< double > &value_in,
compound_type_traits< std::vector< double > >::value_type &value_out)
";

%feature("docstring")  detran_ioutils::get_value "void
detran_ioutils::get_value(T &value_out, const typename
compound_type_traits< T >::value_type &value_in)

Get the value to be stored in memory.

Parameters:
-----------

value_out:  The value to be stored in memory

value_in:  The value as stored in th efile ";

%feature("docstring")  detran_ioutils::get_value "void
detran_ioutils::get_value(std::string &value_out, const
compound_type_traits< std::string >::value_type &value_in) ";

%feature("docstring")  detran_ioutils::get_value "void
detran_ioutils::get_value(std::vector< int > &value_out, const
compound_type_traits< std::vector< int > >::value_type &value_in) ";

%feature("docstring")  detran_ioutils::get_value "void
detran_ioutils::get_value(std::vector< double > &value_out, const
compound_type_traits< std::vector< double > >::value_type &value_in)
";

%feature("docstring")  detran_ioutils::HDF5_MemoryType::type<
std::string > " hid_t detran_ioutils::HDF5_MemoryType::type<
std::string >() ";

%feature("docstring")  detran_ioutils::HDF5_MemoryType::type<
std::vector< int > > " hid_t detran_ioutils::HDF5_MemoryType::type<
std::vector< int > >() ";

%feature("docstring")  detran_ioutils::HDF5_MemoryType::type<
std::vector< double > > " hid_t detran_ioutils::HDF5_MemoryType::type<
std::vector< double > >() ";

%feature("docstring")  detran_ioutils::HDF5_FileType::type<
std::string > " hid_t detran_ioutils::HDF5_FileType::type< std::string
>() ";

%feature("docstring")  detran_ioutils::HDF5_FileType::type<
std::vector< int > > " hid_t detran_ioutils::HDF5_FileType::type<
std::vector< int > >() ";

%feature("docstring")  detran_ioutils::HDF5_FileType::type<
std::vector< double > > " hid_t detran_ioutils::HDF5_FileType::type<
std::vector< double > >() ";

%feature("docstring")  detran_ioutils::print_vec "void
detran_ioutils::print_vec(const T &v, const std::string name=\"\") ";

%feature("docstring")  detran_ioutils::print_vec "void
detran_ioutils::print_vec(const detran_utilities::vec_dbl &v, const
std::string name) ";


// File: namespacedetran__material.xml


// File: namespacedetran__ortho.xml


// File: namespacedetran__postprocess.xml


// File: namespacedetran__test.xml
%feature("docstring")  detran_test::quadruplerange_fixture "static
SP_quadrature detran_test::quadruplerange_fixture()

Create a QuadrupleRange quadrature for testing.

Note, for classes like this, it's nearly as easy simply to create
directly when required. However, setting it once here ensures
consistency at all use points. ";

%feature("docstring")  detran_test::gausslegendre_fixture "static
SP_quadrature detran_test::gausslegendre_fixture()

Create a GaussLegendre quadrature for testing.

Note, for classes like this, it's nearly as easy simply to create
directly when required. However, setting it once here ensures
consistency at all use points. ";

%feature("docstring")  detran_test::uniform_fixture "static
SP_quadrature detran_test::uniform_fixture()

Create a Uniform MOC quadrature for testing.

Note, for classes like this, it's nearly as easy simply to create
directly when required. However, setting it once here ensures
consistency at all use points. ";

%feature("docstring")  detran_test::constant_source_fixture "static
SP_externalsource detran_test::constant_source_fixture(int dimension)
";

%feature("docstring")  detran_test::mesh_1d_fixture "static SP_mesh
detran_test::mesh_1d_fixture()

Create a common 1D mesh for transport tests. ";

%feature("docstring")  detran_test::mesh_2d_fixture "static SP_mesh
detran_test::mesh_2d_fixture(int id=0)

Create a common 2D mesh for transport tests.

This is a square mesh with 2x2 coarse meshes of width 10 cm sides,
each with 10x10 fine meshes. The bottom left coarse mesh has material
0, while the rest have material 1. ";

%feature("docstring")  detran_test::mesh_3d_fixture "static SP_mesh
detran_test::mesh_3d_fixture()

Create a common 3D mesh for transport tests. ";

%feature("docstring")  detran_test::pincell_fixture "static
SP_pincell detran_test::pincell_fixture() ";

%feature("docstring")  detran_test::assembly_fixture "static
SP_assembly detran_test::assembly_fixture() ";

%feature("docstring")  detran_test::core_fixture "static SP_core
detran_test::core_fixture() ";

%feature("docstring")  detran_test::material_fixture_1g "static
SP_material detran_test::material_fixture_1g()

Two group material database for transport tests.

Reference: Mosher's PhD thesis. ";

%feature("docstring")  detran_test::material_fixture_2g "static
SP_material detran_test::material_fixture_2g()

Two group material database for transport tests.

Reference: Mosher's PhD thesis. ";

%feature("docstring")  detran_test::material_fixture_7g "static
SP_material detran_test::material_fixture_7g()

Seven group material database for transport tests.

Reference: C5G7 ";

%feature("docstring")  detran_test::coarsemesh_1d "detran::CoarseMesh::SP_coarsemesh detran_test::coarsemesh_1d() ";

%feature("docstring")  detran_test::coarsemesh_2d "detran::CoarseMesh::SP_coarsemesh detran_test::coarsemesh_2d() ";

%feature("docstring")  detran_test::coarsemesh_2d_b "detran::CoarseMesh::SP_coarsemesh detran_test::coarsemesh_2d_b() ";

%feature("docstring")  detran_test::coarsemesh_3d "detran::CoarseMesh::SP_coarsemesh detran_test::coarsemesh_3d() ";

%feature("docstring")  detran_test::coarsemesh_3d_b "detran::CoarseMesh::SP_coarsemesh detran_test::coarsemesh_3d_b() ";

%feature("docstring")  detran_test::itoa "static std::string
detran_test::itoa(int i)

integer to string ";

%feature("docstring")  detran_test::dtoa "static std::string
detran_test::dtoa(double d)

double to string ";


// File: namespacedetran__user.xml
%feature("docstring")  detran_user::update_T_rhs "void
detran_user::update_T_rhs(void *data, detran::TimeStepper< D > *step,
double t, double dt) ";


// File: namespacedetran__utilities.xml
%feature("docstring")  detran_utilities::itoa "std::string
detran_utilities::itoa(int i) ";

%feature("docstring")  detran_utilities::dtoa "std::string
detran_utilities::dtoa(double d) ";

%feature("docstring")  detran_utilities::InputDB::get< std::string > "
std::string detran_utilities::InputDB::get< std::string >(const
std::string &key) const ";

%feature("docstring")  detran_utilities::InputDB::get<
InputDB::SP_input > " InputDB::SP_input
detran_utilities::InputDB::get< InputDB::SP_input >(const std::string
&key) const ";

%feature("docstring")  detran_utilities::InputDB::get_map< std::string
> " const std::map<std::string, std::string>&
detran_utilities::InputDB::get_map< std::string >() ";

%feature("docstring")  detran_utilities::InputDB::get_map<
InputDB::SP_input > " const std::map<std::string, InputDB::SP_input>&
detran_utilities::InputDB::get_map< InputDB::SP_input >() ";

%feature("docstring")  detran_utilities::norm "double
detran_utilities::norm(vec_dbl &x, std::string flag=\"L2\")

Norm of a vector.

Parameters:
-----------

flag:  L2 (Default), L1, or Linf ";

%feature("docstring")  detran_utilities::vec_scale "void
detran_utilities::vec_scale(vec_dbl &x, double scale)

Scale a double vector. ";

%feature("docstring")  detran_utilities::norm_residual "double
detran_utilities::norm_residual(vec_dbl &x, vec_dbl &y, std::string
flag=\"L2\")

Norm of the residual of two double vectors.

Parameters:
-----------

flag:  Default is false for L2. True for L-infinity. ";

%feature("docstring")  detran_utilities::norm_relative_residual "double detran_utilities::norm_relative_residual(vec_dbl &x, vec_dbl
&y, std::string flag=\"L2\") ";

%feature("docstring")  detran_utilities::distance "double
detran_utilities::distance(Point p1, Point p2)

Distance between two points. ";

%feature("docstring")  detran_utilities::soft_equiv "bool
detran_utilities::soft_equiv(const FPT &value, const FPT &reference,
const FPT precision=1.0e-12)

Compare two floating point scalars for equivalence to a specified
tolerance.

Parameters:
-----------

value:  scalar floating point value

reference:  scalar floating point reference to which value is compared

precision:  tolerance of relative error (default 1.0e-12)

true if values are the same within relative error specified by
precision, false if otherwise Todo Should we be using numeric_limits
instead of hard coded vales for e-12 and e-14? ";

%feature("docstring")  detran_utilities::soft_equiv "bool
detran_utilities::soft_equiv(const int &value, const int &reference,
const int precision) ";

%feature("docstring")  detran_utilities::soft_equiv "bool
detran_utilities::soft_equiv(Value_Iterator value, Value_Iterator
value_end, Ref_Iterator ref, Ref_Iterator ref_end, const typename
std::iterator_traits< Value_Iterator >::value_type precision=1.0e-12)

Compare two floating point fields for equivalence to a specified
tolerance.

Parameters:
-----------

value:  floating point field of values

reference:  floating point field to which values are compared

precision:  tolerance of relative error (default 1.0e-12)

true if values are the same within relative error specified by
precision and the fields are the same size, false if otherwise  The
field soft_equiv check is an element-by-element check of two single-
dimension fields. The precision is the same type as the value field.
The value and reference fields must have STL-type iterators. The
value-types of both fields must be the same or a compile-time error
will result. ";

%feature("docstring")  detran_utilities::warning "void
detran_utilities::warning(int type, std::string message) ";


// File: namespaceexternal__source.xml


// File: namespacegeometry.xml


// File: namespaceioutils.xml


// File: namespacekinetics.xml


// File: namespacematerial.xml


// File: namespaceMGCMTSA.xml


// File: namespaceMGDSA.xml


// File: namespacepostprocess.xml


// File: namespacepython.xml


// File: namespacepython_1_1pydetranutils.xml


// File: namespacepython_1_1pydetranutils_1_1convergence__plot.xml
%feature("docstring")
python::pydetranutils::convergence_plot::plot_eigen_convergence "def
python::pydetranutils::convergence_plot::plot_eigen_convergence Plots
the eigenvalue and fission source error per iteration ";

%feature("docstring")
python::pydetranutils::convergence_plot::plot_convergence "def
python::pydetranutils::convergence_plot::plot_convergence Plots the
residuals from a linear solve ";


// File: namespacepython_1_1pydetranutils_1_1mesh__plot.xml
%feature("docstring")
python::pydetranutils::mesh_plot::plot_mesh_function "def
python::pydetranutils::mesh_plot::plot_mesh_function Plot a mesh
function. ";

%feature("docstring")  python::pydetranutils::mesh_plot::plot_mesh_map
"def python::pydetranutils::mesh_plot::plot_mesh_map Plot a mesh map,
optionally with edges explicitly displayed. ";

%feature("docstring")
python::pydetranutils::mesh_plot::plot_multigroup_flux "def
python::pydetranutils::mesh_plot::plot_multigroup_flux Plot the
multigroup fluxes.        For 1D, they are superimposed on one plot.
In 2D, they     are split into subfigures for the number of groups.
Obviously,     this can get cumbersome for many groups, so we kill it
at 5+. ";

%feature("docstring")  python::pydetranutils::mesh_plot::mesh_axes "def python::pydetranutils::mesh_plot::mesh_axes Get the fine mesh
points for plotting. ";


// File: namespacepython_1_1pydetranutils_1_1quad__plot.xml
%feature("docstring")
python::pydetranutils::quad_plot::plot_quadrature "def
python::pydetranutils::quad_plot::plot_quadrature Plots a quadrature.
";


// File: namespacesolvers.xml


// File: namespacestd.xml


// File: namespacetransport.xml


// File: namespaceutilities.xml


// File: namespaceWGDiffusionLossOperator.xml


// File: ____init_____8py.xml


// File: pydetranutils_2____init_____8py.xml


// File: Acceleration_8cc.xml


// File: Acceleration_8hh.xml


// File: Acceleration_8i_8hh.xml


// File: Assembly_8cc.xml


// File: Assembly_8hh.xml


// File: BDFCoefficients_8hh.xml


// File: BoundaryBase_8hh.xml


// File: BoundaryCondition_8hh.xml


// File: BoundaryConditionMOC_8hh.xml


// File: BoundaryDiffusion_8cc.xml


// File: BoundaryDiffusion_8hh.xml


// File: BoundaryDiffusion_8i_8hh.xml


// File: BoundaryMOC_8cc.xml


// File: BoundaryMOC_8hh.xml


// File: BoundaryMOC_8i_8hh.xml


// File: BoundarySN_8cc.xml


// File: BoundarySN_8hh.xml


// File: BoundarySN_8i_8hh.xml


// File: BoundarySource_8hh.xml


// File: BoundaryTally_8cc.xml


// File: BoundaryTally_8hh.xml


// File: BoundaryTraits_8hh.xml


// File: callow__config_8hh.xml


// File: CallowDefinitions_8hh.xml


// File: ChebyshevDPN_8cc.xml


// File: ChebyshevDPN_8hh.xml


// File: ChebyshevLegendre_8cc.xml


// File: ChebyshevLegendre_8hh.xml


// File: CMR_8cc.xml


// File: CMR_8hh.xml


// File: CMR_8i_8hh.xml


// File: CoarseMesh_8cc.xml


// File: CoarseMesh_8hh.xml


// File: coarsemesh__fixture_8hh.xml


// File: Collocated_8cc.xml


// File: Collocated_8hh.xml


// File: Constants_8hh.xml


// File: ConstantSource_8cc.xml


// File: ConstantSource_8hh.xml


// File: convergence__plot_8py.xml


// File: Core_8cc.xml


// File: Core_8hh.xml


// File: CurrentTally_8cc.xml


// File: CurrentTally_8hh.xml


// File: CurrentTally_8i_8hh.xml


// File: DBC_8hh.xml


// File: Definitions_8hh.xml


// File: detran_8cc.xml
%feature("docstring")  print_welcome "void print_welcome() ";

%feature("docstring")  main "int main(int argc, char **argv) ";


// File: detran__angle_8hh.xml


// File: detran__geometry_8hh.xml


// File: detran__material_8hh.xml


// File: detran__utilities_8hh.xml


// File: DiffusionEigensolver_8cc.xml


// File: DiffusionEigensolver_8hh.xml


// File: DiffusionGainOperator_8cc.xml


// File: DiffusionGainOperator_8hh.xml


// File: DiffusionLossOperator_8cc.xml


// File: DiffusionLossOperator_8hh.xml


// File: DimensionTraits_8hh.xml


// File: DiscreteSource_8cc.xml


// File: DiscreteSource_8hh.xml


// File: DPN_8cc.xml


// File: DPN_8hh.xml


// File: DTN_8cc.xml


// File: DTN_8hh.xml


// File: EigenArnoldi_8cc.xml


// File: EigenArnoldi_8hh.xml


// File: EigenArnoldi_8i_8hh.xml


// File: EigenDiffusion_8cc.xml


// File: EigenDiffusion_8hh.xml


// File: EigenPI_8cc.xml


// File: EigenPI_8hh.xml


// File: EigenPI_8i_8hh.xml


// File: eigenproblem__fixture_8hh.xml


// File: Eigensolver_8cc.xml


// File: EigenSolver_8cc.xml


// File: Eigensolver_8hh.xml


// File: EigenSolver_8hh.xml


// File: EigenSolver_8i_8hh.xml


// File: EigenSolverCreator_8cc.xml


// File: EigenSolverCreator_8hh.xml


// File: EigenvalueManager_8cc.xml


// File: EigenvalueManager_8hh.xml


// File: EnergyIndependentEigenOperator_8cc.xml


// File: EnergyIndependentEigenOperator_8hh.xml


// File: Equation_8hh.xml


// File: Equation__DD__1D_8cc.xml


// File: Equation__DD__1D_8hh.xml


// File: Equation__DD__1D_8i_8hh.xml


// File: Equation__DD__2D_8cc.xml


// File: Equation__DD__2D_8hh.xml


// File: Equation__DD__2D_8i_8hh.xml


// File: Equation__DD__3D_8cc.xml


// File: Equation__DD__3D_8hh.xml


// File: Equation__DD__3D_8i_8hh.xml


// File: Equation__ES__1D_8hh.xml


// File: Equation__MOC_8hh.xml


// File: Equation__SC__1D_8cc.xml


// File: Equation__SC__1D_8hh.xml


// File: Equation__SC__1D_8i_8hh.xml


// File: Equation__SC__2D_8cc.xml


// File: Equation__SC__2D_8hh.xml


// File: Equation__SC__2D_8i_8hh.xml


// File: Equation__SC__MOC_8cc.xml


// File: Equation__SC__MOC_8hh.xml


// File: Equation__SC__MOC_8i_8hh.xml


// File: Equation__SD__1D_8cc.xml


// File: Equation__SD__1D_8hh.xml


// File: Equation__SD__1D_8i_8hh.xml


// File: Equation__SD__2D_8cc.xml


// File: Equation__SD__2D_8hh.xml


// File: Equation__SD__2D_8i_8hh.xml


// File: Execute_8cc.xml


// File: Execute_8hh.xml


// File: external__source__fixture_8hh.xml


// File: ExternalSource_8hh.xml


// File: FissionSource_8cc.xml


// File: FissionSource_8hh.xml


// File: FissionSource_8i_8hh.xml


// File: FixedSourceManager_8cc.xml


// File: FixedSourceManager_8hh.xml


// File: GaussChebyshev_8cc.xml


// File: GaussChebyshev_8hh.xml


// File: GaussLegendre_8cc.xml


// File: GaussLegendre_8hh.xml


// File: GaussSeidel_8cc.xml


// File: GaussSeidel_8hh.xml


// File: GaussSeidel_8i_8hh.xml


// File: GenerateGaussChebyshev_8hh.xml


// File: GenerateGaussLegendre_8hh.xml


// File: GenException_8cc.xml


// File: GenException_8hh.xml


// File: GMRES_8cc.xml


// File: GMRES_8hh.xml


// File: GMRES_8i_8hh.xml


// File: Initialization_8hh.xml
%feature("docstring")  callow_initialize "void callow_initialize(int
argc, char **argv)

Initialize external packages, if enabled. ";

%feature("docstring")  callow_finalize "void callow_finalize()

Finalize external packages, if enabled. ";


// File: InputDB_8cc.xml


// File: InputDB_8hh.xml


// File: InputDB_8i_8hh.xml


// File: IO__HDF5_8cc.xml


// File: IO__HDF5_8hh.xml


// File: IO__HDF5_8t_8hh.xml


// File: IO__HDF5__Traits_8hh.xml


// File: IsotropicSource_8cc.xml


// File: IsotropicSource_8hh.xml


// File: Jacobi_8cc.xml


// File: Jacobi_8hh.xml


// File: Jacobi_8i_8hh.xml


// File: JFNK_8hh.xml


// File: KineticsMaterial_8cc.xml


// File: KineticsMaterial_8hh.xml


// File: KineticsMaterial_8i_8hh.xml


// File: LevelSymmetric_8cc.xml


// File: LevelSymmetric_8hh.xml


// File: LinearExternalSource_8cc.xml


// File: LinearExternalSource_8hh.xml


// File: LinearExternalSource_8i_8hh.xml


// File: LinearMaterial_8cc.xml


// File: LinearMaterial_8hh.xml


// File: LinearSolver_8cc.xml


// File: LinearSolver_8hh.xml


// File: LinearSolver_8i_8hh.xml


// File: LinearSolverCreator_8cc.xml


// File: LinearSolverCreator_8hh.xml


// File: LRA_8cc.xml


// File: LRA_8hh.xml


// File: Manager_8hh.xml


// File: Material_8cc.xml


// File: Material_8hh.xml


// File: Material_8i_8hh.xml


// File: material__fixture_8hh.xml


// File: MathUtilities_8hh.xml


// File: Matrix_8cc.xml


// File: Matrix_8hh.xml


// File: Matrix_8i_8hh.xml


// File: matrix__fixture_8hh.xml


// File: MatrixBase_8hh.xml


// File: MatrixDense_8cc.xml


// File: MatrixDense_8hh.xml


// File: MatrixDense_8i_8hh.xml


// File: MatrixShell_8cc.xml


// File: MatrixShell_8hh.xml


// File: MatrixShell_8i_8hh.xml


// File: matrixshell__fixture_8hh.xml


// File: Mesh_8cc.xml


// File: Mesh_8hh.xml


// File: Mesh_8i_8hh.xml


// File: Mesh1D_8cc.xml


// File: Mesh1D_8hh.xml


// File: Mesh2D_8cc.xml


// File: Mesh2D_8hh.xml


// File: Mesh3D_8cc.xml


// File: Mesh3D_8hh.xml


// File: mesh__fixture_8hh.xml


// File: mesh__plot_8py.xml


// File: MeshMOC_8hh.xml


// File: MGCMTSA_8hh.xml


// File: MGDiffusionSolver_8cc.xml


// File: MGDiffusionSolver_8hh.xml


// File: MGDSA_8cc.xml


// File: MGDSA_8hh.xml


// File: MGPreconditioner_8cc.xml


// File: MGPreconditioner_8hh.xml


// File: MGSolver_8cc.xml


// File: MGSolver_8hh.xml


// File: MGSolver_8i_8hh.xml


// File: MGSolverGMRES_8cc.xml


// File: MGSolverGMRES_8hh.xml


// File: MGSolverGMRES_8i_8hh.xml


// File: MGSolverGS_8cc.xml


// File: MGSolverGS_8hh.xml


// File: MGSolverGS_8i_8hh.xml


// File: MGTransportOperator_8cc.xml


// File: MGTransportOperator_8hh.xml


// File: MGTransportSolver_8cc.xml


// File: MGTransportSolver_8hh.xml


// File: MGTransportSolver_8i_8hh.xml


// File: MomentIndexer_8cc.xml


// File: MomentIndexer_8hh.xml


// File: MomentIndexer_8i_8hh.xml


// File: MomentToDiscrete_8cc.xml


// File: MomentToDiscrete_8hh.xml


// File: MomentToDiscrete_8i_8hh.xml


// File: MultiPhysics_8cc.xml


// File: MultiPhysics_8hh.xml


// File: OrthogonalBasis_8hh.xml


// File: PC__DSA_8cc.xml


// File: PC__DSA_8hh.xml


// File: PC__DSA_8i_8hh.xml


// File: PCIdentity_8hh.xml


// File: PCILU0_8cc.xml


// File: PCILU0_8hh.xml


// File: PCILU0_8i_8hh.xml


// File: PCJacobi_8cc.xml


// File: PCJacobi_8hh.xml


// File: PCJacobi_8i_8hh.xml


// File: PCShell_8cc.xml


// File: PCShell_8hh.xml


// File: PetscSolver_8cc.xml


// File: PetscSolver_8hh.xml


// File: PetscSolver_8i_8hh.xml


// File: PinCell_8cc.xml


// File: PinCell_8hh.xml


// File: Point_8hh.xml


// File: PolarQuadrature_8hh.xml


// File: PowerIteration_8cc.xml


// File: PowerIteration_8hh.xml


// File: PowerIteration_8i_8hh.xml


// File: Preconditioner_8hh.xml


// File: Preconditioner_8i_8hh.xml


// File: PreconditionerBase_8hh.xml


// File: PreconditionerMG_8cc.xml


// File: PreconditionerMG_8hh.xml


// File: PreconditionerMG_8i_8hh.xml


// File: Precursors_8cc.xml


// File: Precursors_8hh.xml


// File: Precursors_8i_8hh.xml


// File: ProductQuadrature_8cc.xml


// File: ProductQuadrature_8hh.xml


// File: ProductQuadrature_8i_8hh.xml


// File: Profiler_8hh.xml


// File: PulsedExternalSource_8cc.xml


// File: PulsedExternalSource_8hh.xml


// File: PulsedExternalSource_8i_8hh.xml


// File: PyExecute_8hh.xml


// File: PyTimeDependentMaterial_8cc.xml


// File: PyTimeDependentMaterial_8hh.xml


// File: quad__plot_8py.xml


// File: Quadrature_8cc.xml


// File: Quadrature_8hh.xml


// File: Quadrature_8i_8hh.xml


// File: quadrature__fixture_8hh.xml


// File: QuadratureFactory_8cc.xml


// File: QuadratureFactory_8hh.xml


// File: QuadratureMOC_8cc.xml


// File: QuadratureMOC_8hh.xml


// File: quadraturemoc__fixture_8hh.xml


// File: QuadrupleRange_8cc.xml


// File: QuadrupleRange_8hh.xml


// File: ReactionRates_8cc.xml


// File: ReactionRates_8hh.xml


// File: Reflective_8hh.xml


// File: Reflective_8i_8hh.xml


// File: ReflectiveMOC_8hh.xml


// File: ReflectiveMOC_8i_8hh.xml


// File: ReflectiveSolver_8cc.xml


// File: ReflectiveSolver_8hh.xml


// File: Richardson_8cc.xml


// File: Richardson_8hh.xml


// File: Richardson_8i_8hh.xml


// File: ScatterSource_8cc.xml


// File: ScatterSource_8hh.xml


// File: ScatterSource_8i_8hh.xml


// File: Segment_8hh.xml


// File: SiloOutput_8cc.xml


// File: SiloOutput_8hh.xml


// File: SlepcSolver_8cc.xml


// File: SlepcSolver_8hh.xml


// File: SlepcSolver_8i_8hh.xml


// File: SoftEquivalence_8hh.xml


// File: Solver_8cc.xml


// File: Solver_8hh.xml


// File: SP_8hh.xml


// File: SP_8i_8hh.xml


// File: SphericalHarmonics_8cc.xml


// File: SphericalHarmonics_8hh.xml


// File: State_8cc.xml


// File: State_8hh.xml


// File: State_8i_8hh.xml


// File: StdOutUtils_8hh.xml


// File: StupidParser_8cc.xml


// File: StupidParser_8hh.xml


// File: Sweeper_8cc.xml


// File: Sweeper_8hh.xml


// File: Sweeper_8t_8hh.xml


// File: Sweeper1D_8cc.xml


// File: Sweeper1D_8hh.xml


// File: Sweeper1D_8i_8hh.xml


// File: Sweeper2D_8cc.xml


// File: Sweeper2D_8hh.xml


// File: Sweeper2D_8i_8hh.xml


// File: Sweeper2DMOC_8cc.xml


// File: Sweeper2DMOC_8hh.xml


// File: Sweeper2DMOC_8i_8hh.xml


// File: Sweeper3D_8cc.xml


// File: Sweeper3D_8hh.xml


// File: Sweeper3D_8i_8hh.xml


// File: SweepOperator_8cc.xml


// File: SweepOperator_8hh.xml


// File: SweepSource_8hh.xml


// File: SweepSource_8i_8hh.xml


// File: SyntheticDiscreteSource_8cc.xml


// File: SyntheticDiscreteSource_8hh.xml


// File: SyntheticDiscreteSource_8i_8hh.xml


// File: SyntheticMomentSource_8cc.xml


// File: SyntheticMomentSource_8hh.xml


// File: SyntheticMomentSource_8i_8hh.xml


// File: SyntheticSource_8cc.xml


// File: SyntheticSource_8hh.xml


// File: SyntheticSource_8i_8hh.xml


// File: TabuchiYamamoto_8cc.xml


// File: TabuchiYamamoto_8hh.xml


// File: test__BoundaryDiffusion_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_BoundaryDiffusion_1D "int
test_BoundaryDiffusion_1D(int argc, char *argv[]) ";

%feature("docstring")  test_BoundaryDiffusion_2D "int
test_BoundaryDiffusion_2D(int argc, char *argv[]) ";

%feature("docstring")  test_BoundaryDiffusion_3D "int
test_BoundaryDiffusion_3D(int argc, char *argv[]) ";


// File: test__BoundaryMOC_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_BoundaryMOC "int test_BoundaryMOC(int
argc, char *argv[]) ";

%feature("docstring")  test_BoundaryMOC_2 "int test_BoundaryMOC_2(int
argc, char *argv[]) ";


// File: test__BoundarySN_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_BoundarySN "int test_BoundarySN(int argc,
char *argv[]) ";


// File: test__CoarseMesh_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_CoarseMesh "int test_CoarseMesh(int argc,
char *argv[]) ";


// File: test__Collocated_8cc.xml
%feature("docstring")  std::main "int main(int argc, char *argv[]) ";

%feature("docstring")  std::test_Collocated "int test_Collocated(int
argc, char *argv[]) ";


// File: test__ConstantSource_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_ConstantSource "int
test_ConstantSource(int argc, char *argv[]) ";


// File: test__CurrentTally_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_CurrentTally_1D "int
test_CurrentTally_1D(int argc, char *argv[]) ";

%feature("docstring")  test_CurrentTally_2D "int
test_CurrentTally_2D(int argc, char *argv[]) ";

%feature("docstring")  test_CurrentTally_3D "int
test_CurrentTally_3D(int argc, char *argv[]) ";


// File: test__DiffusionEigensolver_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_DiffusionEigensolver "int
test_DiffusionEigensolver(int argc, char *argv[]) ";


// File: test__DiffusionFixedSourceSolver_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_DiffusionFixedSourceSolver_1D "int
test_DiffusionFixedSourceSolver_1D(int argc, char *argv[]) ";


// File: test__DiscreteSource_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_DiscreteSource "int
test_DiscreteSource(int argc, char *argv[]) ";


// File: test__EigenSolver_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_PowerIteration "int
test_PowerIteration(int argc, char *argv[]) ";

%feature("docstring")  test_SlepcSolver "int test_SlepcSolver(int
argc, char *argv[]) ";


// File: test__EigenvalueManager_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_EigenvalueManager_material "SP_material
test_EigenvalueManager_material() ";

%feature("docstring")  test_EigenvalueManager_mesh "SP_mesh
test_EigenvalueManager_mesh(int d) ";

%feature("docstring")  test_EigenvalueManager_input "InputDB::SP_input test_EigenvalueManager_input() ";

%feature("docstring")  test_EigenvalueManager_T "int
test_EigenvalueManager_T() ";

%feature("docstring")  test_EigenvalueManager_1D "int
test_EigenvalueManager_1D(int argc, char *argv[]) ";

%feature("docstring")  test_EigenvalueManager_2D "int
test_EigenvalueManager_2D(int argc, char *argv[]) ";

%feature("docstring")  test_EigenvalueManager_3D "int
test_EigenvalueManager_3D(int argc, char *argv[]) ";


// File: test__Equation__DD__1D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Equation_DD_1D "int
test_Equation_DD_1D(int argc, char *argv[]) ";


// File: test__Equation__DD__2D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Equation_DD_2D_basic "int
test_Equation_DD_2D_basic(int argc, char *argv[]) ";


// File: test__Equation__SC__1D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Equation_SC_1D "int
test_Equation_SC_1D(int argc, char *argv[]) ";


// File: test__Equation__SC__MOC_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Equation_SC_MOC "int
test_Equation_SC_MOC(int argc, char *argv[]) ";


// File: test__FixedSourceManager_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_FixedSourceManager_material "SP_material
test_FixedSourceManager_material() ";

%feature("docstring")  test_FixedSourceManager_mesh "SP_mesh
test_FixedSourceManager_mesh(int d) ";

%feature("docstring")  test_FixedSourceManager_input "InputDB::SP_input test_FixedSourceManager_input() ";

%feature("docstring")  compute_leakage "double
compute_leakage(typename BoundarySN< D >::SP_boundary b, SP_quadrature
q, SP_material mat, SP_mesh mesh) ";

%feature("docstring")  test_FixedSourceManager_T "int
test_FixedSourceManager_T() ";

%feature("docstring")  test_FixedSourceManager_1D "int
test_FixedSourceManager_1D(int argc, char *argv[]) ";

%feature("docstring")  test_FixedSourceManager_2D "int
test_FixedSourceManager_2D(int argc, char *argv[]) ";

%feature("docstring")  test_FixedSourceManager_3D "int
test_FixedSourceManager_3D(int argc, char *argv[]) ";

%feature("docstring")  test_FixedSourceManager_iterate "int
test_FixedSourceManager_iterate(int argc, char *argv[]) ";


// File: test__GaussLegendre_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_GaussLegendre_basic "int
test_GaussLegendre_basic(int argc, char *argv[]) ";


// File: test__InnerGMRES_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_InnerGMRES_1D_actual "int
test_InnerGMRES_1D_actual() ";

%feature("docstring")  test_InnerGMRES_1D "int test_InnerGMRES_1D(int
argc, char *argv[]) ";


// File: test__InputDB_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_InputDB "int test_InputDB(int argc, char
*argv[]) ";

%feature("docstring")  test_InputDB_serialize "int
test_InputDB_serialize(int argc, char *argv[]) ";


// File: test__IO__HDF5_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_IO_HDF5_input "int test_IO_HDF5_input(int
argc, char *argv[]) ";

%feature("docstring")  test_IO_HDF5_material "int
test_IO_HDF5_material(int argc, char *argv[]) ";

%feature("docstring")  test_IO_HDF5_mesh "int test_IO_HDF5_mesh(int
argc, char *argv[]) ";


// File: test__IsotropicSource_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_IsotropicSource "int
test_IsotropicSource(int argc, char *argv[]) ";


// File: test__KineticsMaterial_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_KineticsMaterial "int
test_KineticsMaterial(int argc, char *argv[]) ";


// File: test__LinearExternalSource_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_LinearExternalSource "int
test_LinearExternalSource(int argc, char *argv[]) ";


// File: test__LinearMaterial_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_LinearMaterial "int
test_LinearMaterial(int argc, char *argv[]) ";


// File: test__LinearSolver_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  get_db "LinearSolver::SP_db get_db() ";

%feature("docstring")  test_Richardson "int test_Richardson(int argc,
char *argv[]) ";

%feature("docstring")  test_Jacobi "int test_Jacobi(int argc, char
*argv[]) ";

%feature("docstring")  test_GaussSeidel "int test_GaussSeidel(int
argc, char *argv[]) ";

%feature("docstring")  test_SOR "int test_SOR(int argc, char *argv[])
";

%feature("docstring")  test_GMRES "int test_GMRES(int argc, char
*argv[]) ";

%feature("docstring")  test_PetscSolver "int test_PetscSolver(int
argc, char *argv[]) ";


// File: test__LRA_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_monitor "void test_monitor(void *data,
TimeStepper< _2D > *ts, int step, double t, double dt, int it, bool
conv) ";

%feature("docstring")  get_mesh "Mesh2D::SP_mesh
get_mesh(Mesh2D::size_t fmm=1) ";

%feature("docstring")  test_LRA "int test_LRA(int argc, char *argv[])
";


// File: test__Material_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Material_basic "int
test_Material_basic(int argc, char *argv[]) ";

%feature("docstring")  test_Material_bounds "int
test_Material_bounds(int argc, char *argv[]) ";

%feature("docstring")  test_Material_serialize "int
test_Material_serialize(int argc, char *argv[]) ";


// File: test__MathUtilities_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_norm_L2 "int test_norm_L2(int argc, char
*argv[]) ";

%feature("docstring")  test_norm_L1 "int test_norm_L1(int argc, char
*argv[]) ";

%feature("docstring")  test_norm_Linf "int test_norm_Linf(int argc,
char *argv[]) ";

%feature("docstring")  test_vec_scale "int test_vec_scale(int argc,
char *argv[]) ";

%feature("docstring")  test_norm_residual_L2 "int
test_norm_residual_L2(int argc, char *argv[]) ";

%feature("docstring")  test_norm_residual_L1 "int
test_norm_residual_L1(int argc, char *argv[]) ";

%feature("docstring")  test_norm_residual_Linf "int
test_norm_residual_Linf(int argc, char *argv[]) ";


// File: test__Matrix_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Matrix "int test_Matrix(int argc, char
*argv[])

1 2 3 4 5 6 ";


// File: test__MatrixDense_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MatrixDense "int test_MatrixDense(int
argc, char *argv[]) ";


// File: test__MatrixShell_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MatrixShell "int test_MatrixShell(int
argc, char *argv[]) ";


// File: test__Mesh1D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Mesh1D "int test_Mesh1D(int argc, char
*argv[]) ";

%feature("docstring")  test_Mesh1D_serialize "int
test_Mesh1D_serialize(int argc, char *argv[]) ";


// File: test__Mesh2D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Mesh2D_basic "int test_Mesh2D_basic(int
argc, char *argv[]) ";


// File: test__Mesh3D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Mesh3D "int test_Mesh3D(int argc, char
*argv[]) ";


// File: test__MomentIndexer_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MomentIndexer "int test_MomentIndexer(int
argc, char *argv[]) ";


// File: test__MomentToDiscrete_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MomentToDiscrete "int
test_MomentToDiscrete(int argc, char *argv[]) ";


// File: test__MultiPhysics_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MultiPhysics "int test_MultiPhysics(int
argc, char *argv[]) ";


// File: test__PinCell_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_PinCell "int test_PinCell(int argc, char
*argv[]) ";


// File: test__Point_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Point "int test_Point(int argc, char
*argv[]) ";


// File: test__PowerIteration_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_PowerIteration_2D "int
test_PowerIteration_2D(int argc, char *argv[]) ";


// File: test__PulsedExternalSource_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_PulsedExternalSource "int
test_PulsedExternalSource(int argc, char *argv[]) ";


// File: test__QuadrupleRange_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_QuadrupleRange_basic "int
test_QuadrupleRange_basic(int argc, char *argv[]) ";


// File: test__ReactionRates_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_ReactionRates_pinpower "int
test_ReactionRates_pinpower(int argc, char *argv[]) ";


// File: test__Segment_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Segment "int test_Segment(int argc, char
*argv[]) ";


// File: test__SiloOutput_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_SiloOutput "int test_SiloOutput(int argc,
char *argv[]) ";


// File: test__SourceIteration_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_SourceIteration_2D "int
test_SourceIteration_2D(int argc, char *argv[]) ";


// File: test__SphericalHarmonics_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_SphericalHarmonics "int
test_SphericalHarmonics(int argc, char *argv[]) ";

%feature("docstring")  test_SphericalHarmonics_integration "int
test_SphericalHarmonics_integration(int argc, char *argv[]) ";


// File: test__State_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_State_basic "int test_State_basic(int
argc, char *argv[]) ";


// File: test__StupidParser_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_StupidParser "int test_StupidParser(int
argc, char *tmp_argv[]) ";

%feature("docstring")  test_StupidParser_hdf5 "int
test_StupidParser_hdf5(int argc, char *tmp_argv[]) ";


// File: test__Sweeper2D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Sweeper2D_basic "int
test_Sweeper2D_basic(int argc, char *argv[]) ";


// File: test__Sweeper3D_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Sweeper3D_basic "int
test_Sweeper3D_basic(int argc, char *argv[]) ";


// File: test__TabuchiYamamoto_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_TabuchiYamamoto "int
test_TabuchiYamamoto(int argc, char *argv[]) ";


// File: test__Testing_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Testing_pass "int test_Testing_pass(int
argc, char *argv[]) ";

%feature("docstring")  test_Testing_fail "int test_Testing_fail(int
argc, char *argv[]) ";


// File: test__TimeStepper_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_monitor "void test_monitor(void *data,
TimeStepper< _1D > *ts, int step, double t, double dt, int it, bool
conv) ";

%feature("docstring")  test_TimeStepper "int test_TimeStepper(int
argc, char *argv[])

1D time dependent SN problem with a TD source and no fission. Constant
source from 0:1. Final time 10. ";

%feature("docstring")  test_BDF_Steps "int test_BDF_Steps(int argc,
char *argv[])

Test BDF steps agains reference. ";


// File: test__Track_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Track "int test_Track(int argc, char
*argv[]) ";


// File: test__Tracker_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Tracker_2x2 "int test_Tracker_2x2(int
argc, char *argv[]) ";

%feature("docstring")  test_Tracker_3x3 "int test_Tracker_3x3(int
argc, char *argv[]) ";


// File: test__TWIGL_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_monitor "void test_monitor(void *data,
TimeStepper< _2D > *ts, int step, double t, double dt, int it, bool
conv) ";

%feature("docstring")  get_mesh "Mesh2D::SP_mesh
get_mesh(Mesh2D::size_t fmm=1) ";

%feature("docstring")  test_TWIGL "int test_TWIGL(int argc, char
*argv[]) ";


// File: test__Uniform_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Uniform "int test_Uniform(int argc, char
*argv[]) ";


// File: test__Vector_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Vector "int test_Vector(int argc, char
*argv[]) ";

%feature("docstring")  test_Vector_resize "int test_Vector_resize(int
argc, char *argv[]) ";


// File: TestDriver_8hh.xml


// File: TimeDependentExternalSource_8hh.xml


// File: TimeDependentManager_8cc.xml


// File: TimeDependentManager_8hh.xml


// File: TimeDependentMaterial_8cc.xml


// File: TimeDependentMaterial_8hh.xml


// File: Timer_8hh.xml


// File: TimeStepper_8cc.xml


// File: TimeStepper_8hh.xml


// File: Track_8hh.xml


// File: TrackDB_8cc.xml


// File: TrackDB_8hh.xml


// File: Tracker_8cc.xml


// File: Tracker_8hh.xml


// File: TransportManager_8hh.xml


// File: Triplet_8hh.xml


// File: Typedefs_8hh.xml


// File: Uniform_8cc.xml


// File: Uniform_8hh.xml


// File: UniformEqual_8cc.xml


// File: UniformEqual_8hh.xml


// File: Vacuum_8hh.xml


// File: VacuumMOC_8hh.xml


// File: Vector_8cc.xml


// File: Vector_8hh.xml


// File: Vector_8i_8hh.xml


// File: Warning_8hh.xml


// File: WGDiffusionLossOperator_8cc.xml


// File: WGDiffusionLossOperator_8hh.xml


// File: WGPreconditioner_8cc.xml


// File: WGPreconditioner_8hh.xml


// File: WGSolver_8cc.xml


// File: WGSolver_8hh.xml


// File: WGSolver_8t_8hh.xml


// File: WGSolverGMRES_8cc.xml


// File: WGSolverGMRES_8hh.xml


// File: WGSolverGMRES_8i_8hh.xml


// File: WGSolverSI_8cc.xml


// File: WGSolverSI_8hh.xml


// File: WGSolverSI_8i_8hh.xml


// File: WGTransportOperator_8cc.xml


// File: WGTransportOperator_8hh.xml


// File: WithinGroupAcceleration_8hh.xml


// File: callow.xml


// File: DBC.xml


// File: testing.xml


// File: todo.xml


// File: dir_ea014b4352ba9f6df6b15a0d8085a64d.xml


// File: dir_12ca6eed5f6d0ea8ad0b8c8baf252f9e.xml


// File: dir_0c18ac3aba294f288685bdffe60fef72.xml


// File: dir_d836f0bd14e2f99a98c797f06920e049.xml


// File: dir_0202c16cfac47ada0ad714fa9c6e3b04.xml


// File: dir_35f9d6256392071ebc3bd33bf80d6afd.xml


// File: dir_b879090c8e9cd8d9e8a78b3d8600c7d2.xml


// File: dir_9355da5016324d2c3af067cb3872c636.xml


// File: dir_b5a48f68a85ed86b439f40a9cfc7818d.xml


// File: dir_cdb218c2bca9f1924d7f10b4b803813e.xml


// File: dir_858f8440659bb7580eb4a23721f77640.xml


// File: dir_9bbdcfdba5aef14129ee854316e38c8e.xml


// File: dir_1e115a3d91e459290bdcd2790be1d362.xml


// File: dir_28734ae0835b5d8e9fa9bd0c2998effd.xml


// File: dir_5ce853c0c2c6b66ec969619743521c20.xml


// File: dir_5ebd9147cf49bc4c69e95e60a8dcef9e.xml


// File: dir_33d5f474ec21ac5009ab645e3ae44662.xml


// File: dir_b62eefc9085a412cc1c8f53f017bb9ea.xml


// File: dir_44f7a74b83fe8f09fd24d2a1bc2cfd93.xml


// File: dir_2f512e57ee91b3828f92b485b05e9413.xml


// File: dir_eeb5eb01ba66297d3763a4e3b5d0f060.xml


// File: dir_2f78a55c3cba4bb969556d5ad7872390.xml


// File: dir_f1fd10f168e4359a664b19a2ac82329c.xml


// File: dir_0266203dbeea4ff8942a20e0b10d485f.xml


// File: dir_f9f7e1d329c10a6c7b3a426a3f50ef04.xml


// File: dir_0f2313b097fccbc913ea9f1f0ffc45e3.xml


// File: dir_c0c950a5f879b845fc41b15d360fbe18.xml


// File: dir_33e7fb1b1271bd1a4d46130da11d2edf.xml


// File: dir_3fc1e4be4484d4cd0bed8658c904e9e8.xml


// File: dir_9bc16325138b74d5895f17576ea57b40.xml


// File: dir_763f5802e192b12f56a80516023abcb2.xml


// File: dir_0d29a023a4197274ab6b1884361d5782.xml


// File: dir_454554df072bd86a548228f0628ccb39.xml


// File: dir_0f87e35d544d46f3aea4e3adeecb0fab.xml


// File: dir_9e2a030e98e845d5693be766170ee212.xml


// File: dir_aed5e7f73191ff854208fb501dbd009b.xml


// File: dir_a3091ec8ae2a2d761a97dfebd2b41cb8.xml


// File: dir_79d5b10765767b35747a6f59fdb9c276.xml


// File: dir_721682c634a8df44382245b1e97edd7c.xml


// File: _2transport_2test_2test_CurrentTally_8cc-example.xml


// File: angle_2test_2test_QuadrupleRange_8cc-example.xml


// File: angle_2test_2test_SphericalHarmonics_8cc-example.xml


// File: geometry_2test_2test_TrackDB_8cc-example.xml


// File: test_2test_MomentToDiscrete_8cc-example.xml


// File: utilities_2test_2test_SP_8cc-example.xml


// File: indexpage.xml

