//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   FixedSourceManager.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  FixedSourceManager class definition.
 */
//---------------------------------------------------------------------------//

#include "FixedSourceManager.hh"
#include "SolverMG.hh"
#include "config/detran_config.hh"
#include "angle/QuadratureFactory.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundarySN.hh"
#include "geometry/Tracker.hh"

// Multigroup solvers
#include "GaussSeidelMG.hh"
#ifdef DETRAN_ENABLE_PETSC
#include "petsc.h"
#include "KrylovMG.hh"
#endif


#include <string>

namespace detran
{

template <class D>
FixedSourceManager<D>::FixedSourceManager(SP_input    input,
                                          SP_material material,
                                          SP_mesh     mesh)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_fixed_type(FIXED)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Postconditions
}

template <class D>
void FixedSourceManager<D>::setup()
{

  // Decide the discretization.  Default is discrete ordinates
  // with a diamond difference approximation.
  std::string eq = "dd";
  d_discretization = SN;
  if (d_input->check("equation"))
    eq = d_input->get<std::string>("equation");
  if (eq == "scmoc" || eq == "ddmoc")
    d_discretization = MOC;
  else if (eq == "diffusion")
    d_discretization = DIFF;

  // Setup the quadrature if needed
  if (d_discretization != DIFF)
  {
    detran_angle::QuadratureFactory quad_factory;
    quad_factory.build(d_quadrature, d_input, D::dimension);
    Assert(d_quadrature);
    if (d_discretization == MOC)
    {
      // Track the mesh
      detran_geometry::Tracker tracker(d_mesh, d_quadrature);
      // Normalize segments to conserve volume.
      tracker.normalize();
      // Replace the mesh with the tracked one.  This suggests refactoring
      // to have a (possibly null) trackdb in Mesh.
      d_mesh = tracker.meshmoc();
    }
  }

  // Setup the boundary conditions
  if (d_discretization == MOC)
    d_boundary = new BoundaryMOC<D>(d_input, d_mesh, d_quadrature);
  else if (d_discretization == SN)
    d_boundary = new BoundarySN<D>(d_input, d_mesh, d_quadrature);
  else
    d_boundary = new BoundaryDiffusion<D>(d_input, d_mesh);

  // Setup the state vector.
  // \todo State needs only one constructor with optional quadrature arg
  if (d_discretization == DIFF)
    d_state = new State(d_input, d_mesh);
  else
    d_state = new State(d_input, d_mesh, d_quadrature);

  // Build the fission source if needed
  if (d_input->check("problem_type"))
  {
    if(d_input->get<std::string>("problem_type") == "fixed_multiply")
    {
      d_fissionsource = new FissionSource(d_state, d_mesh, d_material);
      d_fixed_type = FIXEDMULTIPLY;
    }
  }

  // Postconditions
  Ensure(d_boundary);
  Ensure(d_state);
}

template <class D>
void FixedSourceManager<D>::set_source(SP_source q)
{
  if (!q)
  {
    std::cout << "The source q is null.  Ignoring it." << std::endl;
    return;
  }
  // Otherwise, add it to the source list.
  d_sources.push_back(q);
}

template <class D>
void FixedSourceManager<D>::solve()
{
  if (d_discretization == DIFF)
  {
    THROW("NOT FOR DIFFUSION YET");
  }
  else
  {
    typename SolverMG<D>::SP_solver solver;

    // Default solver is Gauss-Seidel
    std::string outer_solver = "GS";
    if (d_input->check("outer_solver"))
    {
      outer_solver = d_input->get<std::string>("outer_solver");
    }

    // \todo enable use of more than one source
    if (outer_solver == "GS")
    {
      solver = new GaussSeidelMG<D>(d_input, d_state, d_mesh, d_material, d_quadrature,
                                    d_boundary, d_sources[0], d_fissionsource);
    }
    else if (outer_solver == "KrylovMG")
    {
      #ifdef DETRAN_ENABLE_PETSC
        solver = new KrylovMG<D>(d_input, d_state, d_mesh, d_material, d_quadrature,
                                 boundary, d_sources[0], d_fissionsource);
      #else
        THROW("KrylovMG is unavailable since PETSc is not enabled.");
      #endif
    }
    else
    {
      THROW("Unsupported outer_solver type selected: " + outer_solver);
    }

    // Solve the problem
    solver->solve();
  }

}

} // end namespace detran


