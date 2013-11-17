//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  FixedSourceManager.cc
 *  @brief FixedSourceManager class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "FixedSourceManager.hh"
#include "angle/QuadratureFactory.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryFactory.t.hh"
#include "geometry/Tracker.hh"
// Multigroup solvers
#include "MGSolverGS.hh"
#include "MGDiffusionSolver.hh"
#include "MGSolverGMRES.hh"
#include "MGSolverCMFD.hh"
#include <string>


namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
FixedSourceManager<D>::FixedSourceManager(int          argc,
                                          char        *argv[],
                                          SP_input     input,
                                          SP_material  material,
                                          SP_mesh      mesh,
                                          bool         multiply,
                                          bool         fission)
  : TransportManager(argc, argv)
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_adjoint(false)
  , d_discretization(0)
  , d_multiply(multiply)
  , d_fission(multiply || fission)
  , d_is_setup(false)
  , d_is_ready(false)
  , d_generation(0)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
}

//----------------------------------------------------------------------------//
template <class D>
FixedSourceManager<D>::FixedSourceManager(SP_input    input,
                                          SP_material material,
                                          SP_mesh     mesh,
                                          bool        multiply,
                                          bool        fission)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_adjoint(false)
  , d_discretization(0)
  , d_multiply(multiply)
  , d_fission(multiply || fission)
  , d_is_setup(false)
  , d_is_ready(false)
  , d_generation(0)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Insist(d_mesh->dimension() == D::dimension,
         "Solver dimension != and mesh dimension.");
}


//----------------------------------------------------------------------------//
template <class D>
void FixedSourceManager<D>::setup()
{
  //--------------------------------------------------------------------------//
  // DISCRETIZATION
  //--------------------------------------------------------------------------//

  std::string eq = "dd";
  d_discretization = SN;
  if (d_input->check("equation"))
    eq = d_input->get<std::string>("equation");
  if (eq == "scmoc" || eq == "ddmoc")
    d_discretization = MOC;
  else if (eq == "diffusion")
    d_discretization = DIFF;

  //--------------------------------------------------------------------------//
  // QUADRATURE
  //--------------------------------------------------------------------------//

  // Setup the quadrature if needed
  if (d_discretization != DIFF)
  {
    d_quadrature =
      detran_angle::QuadratureFactory::build(d_input, D::dimension);
    Assert(d_quadrature);
    if (d_discretization == MOC)
    {
      // Track the mesh
      detran_geometry::Tracker tracker(d_mesh, d_quadrature);
      // Normalize segments to conserve volume.
      tracker.normalize();
      // Replace the mesh with the tracked one.  This suggests refactoring
      // to have a (possibly null) trackdb in Mesh.
      //d_mesh = tracker.meshmoc();
    }
  }

  //--------------------------------------------------------------------------//
  // BOUNDARY AND STATE
  //--------------------------------------------------------------------------//

  // Setup the boundary conditions
  if (d_discretization == MOC)
  {
    d_boundary = BoundaryFactory<D, BoundaryMOC>::
      build(d_input, d_mesh, d_quadrature);
  }
  else if (d_discretization == SN)
  {
    d_boundary = BoundaryFactory<D, BoundarySN>::
      build(d_input, d_mesh, d_quadrature);
  }
  else
  {
    d_boundary = BoundaryFactory<D, BoundaryDiffusion>::
      build(d_input, d_mesh, d_quadrature);
  }
  Insist(d_boundary, "bad boundary");
  // Setup the state vector.
  d_state = new State(d_input, d_mesh, d_quadrature);

  //--------------------------------------------------------------------------//
  // FISSION SOURCE
  //--------------------------------------------------------------------------//

  if (d_fission)
    d_fissionsource = new FissionSource(d_state, d_mesh, d_material);

  // Signify the manager is ready to solve
  d_is_setup = true;

  // Ensure a solver is rebuilt
  d_is_ready = false;

  Ensure(d_boundary);
  Ensure(d_state);
}

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
template <class D>
bool FixedSourceManager<D>::set_solver()
{
  if (!d_is_setup)
  {
    std::cout << "You must setup the manager before solving.  Skipping solve."
              << std::endl;
    return false;
  }

  if (d_discretization == DIFF)
  {
    // Only option is explicit matrix construction
    d_solver = new MGDiffusionSolver<D>(d_state, d_material, d_boundary,
                                        d_sources, d_fissionsource,
                                        d_multiply);
    d_is_ready = true;
  }
  else
  {
    // Default solver is Gauss-Seidel
    std::string outer_solver = "GS";
    if (d_input->check("outer_solver"))
      outer_solver = d_input->get<std::string>("outer_solver");

    if (outer_solver == "GS")
    {
      d_solver = new MGSolverGS<D>(d_state, d_material, d_boundary,
                                   d_sources, d_fissionsource, d_multiply);
    }
    else if (outer_solver == "CMFD")
    {
      d_solver = new MGSolverCMFD<D>(d_state, d_material, d_boundary,
                                     d_sources, d_fissionsource, d_multiply);
    }
    else if (outer_solver == "GMRES")
    {
      d_solver = new MGSolverGMRES<D>(d_state, d_material, d_boundary,
                                      d_sources, d_fissionsource, d_multiply);
    }
    else
    {
      std::cout << "Unsupported outer_solver type selected:"
                << outer_solver << std::endl;
      return false;
    }
    d_is_ready = true;
  }
  return d_is_ready;
}

//----------------------------------------------------------------------------//
template <class D>
bool FixedSourceManager<D>::solve(const double keff)
{

  if (!d_is_ready)
  {
    std::cout << "You must set the solver before solving.  Skipping solve."
              << std::endl;
    return false;
  }

  // Solve the problem
  d_solver->solve(keff);

  return true;
}

//----------------------------------------------------------------------------//
template <class D>
double FixedSourceManager<D>::iterate(const int generation)
{
  // Preconditions
  Require(!d_multiply);
  Require(d_fissionsource);
  Require(d_is_ready);

  using detran_utilities::norm_residual;

  //  setup problem
  //  build fixed source
  //  solve generation 0
  //  save state & boundary for response
  //  for generation = 1, N
  //    compute first fission source
  //    wipe fixed source & boundary
  //    solve generation
  //    save state & boundary

  // Save current group flux.
  State::group_moments_type phi_old = d_state->all_phi();

  if (generation == 0)
  {
    // Solve the direct problem.
    d_solver->solve();
    // Kill the fixed source vector and rebuild the solver with no
    // external sources.
    d_sources.clear();
    set_solver();
  }
  else
  {
    // Build fission source from current state.
    d_fissionsource->update();
    // Wipe state and boundary
    d_state->clear();
    d_boundary->clear();
    // Solve the problem
    d_solver->solve();
  }

  // Compute L2 norm of difference between iterates
  double norm_delta_phi   = 0;
  double norm_delta_phi_g = 0;
  for (int g = 0; g < d_material->number_groups(); ++g)
  {
    norm_delta_phi_g = norm_residual(d_state->phi(g), phi_old[g], "L2");
    norm_delta_phi += norm_delta_phi_g * norm_delta_phi_g;
  }
  norm_delta_phi = std::sqrt(norm_delta_phi);

  return norm_delta_phi;
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

SOLVERS_INSTANTIATE_EXPORT(FixedSourceManager<_1D>)
SOLVERS_INSTANTIATE_EXPORT(FixedSourceManager<_2D>)
SOLVERS_INSTANTIATE_EXPORT(FixedSourceManager<_3D>)
SOLVERS_TEMPLATE_EXPORT(detran_utilities::SP<FixedSourceManager<_1D> >)
SOLVERS_TEMPLATE_EXPORT(detran_utilities::SP<FixedSourceManager<_2D> >)
SOLVERS_TEMPLATE_EXPORT(detran_utilities::SP<FixedSourceManager<_3D> >)

} // end namespace detran


