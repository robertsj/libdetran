//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  EigenvalueManager.cc
 *  @brief EigenvalueManager class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "EigenvalueManager.hh"
#include "angle/QuadratureFactory.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundarySN.hh"
#include "geometry/Tracker.hh"

// Eigenvalue solvers
#include "solvers/eigen/EigenPI.hh"
#include "solvers/eigen/EigenDiffusion.hh"
#include "solvers/eigen/EigenArnoldi.hh"
#include "solvers/eigen/EigenGD.hh"
#include "solvers/eigen/EigenCMFD.hh"

#include <string>

namespace detran
{

using std::cout;
using std::endl;

//---------------------------------------------------------------------------//
template <class D>
EigenvalueManager<D>::EigenvalueManager(int          argc,
                                        char        *argv[],
                                        SP_input     input,
                                        SP_material  material,
                                        SP_mesh      mesh)
  : TransportManager(argc, argv)
  , d_adjoint(false)
  , d_discretization(0)
  , d_is_setup(false)
{
  Require(input);
  Require(material);
  Require(mesh);

  /// Create the fixed source manager
  d_mg_solver = new FixedSourceManager<D>(input, material, mesh, false, true);
  d_mg_solver->setup();
  d_mg_solver->set_solver();
  d_discretization = d_mg_solver->discretization();

  Ensure(d_mg_solver);
}

//---------------------------------------------------------------------------//
template <class D>
EigenvalueManager<D>::EigenvalueManager(SP_input    input,
                                        SP_material material,
                                        SP_mesh     mesh)
  : d_adjoint(false)
  , d_discretization(0)
  , d_is_setup(false)
{
  Require(input);
  Require(material);
  Require(mesh);

  /// Create the fixed source manager
  d_mg_solver = new FixedSourceManager<D>(input, material, mesh, false, true);
  d_mg_solver->setup();
  d_mg_solver->set_solver();
  d_discretization = d_mg_solver->discretization();

  Ensure(d_mg_solver);
}

//---------------------------------------------------------------------------//
template <class D>
bool EigenvalueManager<D>::solve()
{
  std::cout << "Solving eigenvalue problem..." << std::endl;

  std::string eigen_solver = "PI";
  if (d_mg_solver->input()->check("eigen_solver"))
  {
    eigen_solver =
      d_mg_solver->input()->template get<std::string>("eigen_solver");
  }
  if (eigen_solver == "PI")
  {
    d_solver = new EigenPI<D>(d_mg_solver);
  }
  else if (eigen_solver == "diffusion")
  {
    if (d_discretization != Fixed_T::DIFF)
    {
      std::cout << "Diffusion eigensolver requires diffusion discretization."
                << std::endl;
      return false;
    }
    d_solver = new EigenDiffusion<D>(d_mg_solver);
  }
  else if (eigen_solver == "arnoldi")
  {
      d_solver = new EigenArnoldi<D>(d_mg_solver);
  }
  else if (eigen_solver == "cmfd")
  {
      d_solver = new EigenCMFD<D>(d_mg_solver);
  }
  else if (eigen_solver == "GD")
  {
    std::cout << " GD-----> " << std::endl;
    if (d_discretization == Fixed_T::DIFF)
    {
      cout << "GD not applicable for diffusion.  Use the diffusion" << endl;
      cout << "eigensolver and select gd from callow to use the " << endl;
      cout << "built-in implementation or use the SLEPc version. " << endl;
      return false;
    }
    d_solver = new EigenGD<D>(d_mg_solver);
  }
  else
  {
    std::cout << "Unsupported outer_solver type selected:"
              << eigen_solver << std::endl;
    return false;
  }

  // Solve the eigenvalue problem
  d_solver->solve();
  return true;
}

//---------------------------------------------------------------------------//
// Explicit instantiations
//---------------------------------------------------------------------------//

SOLVERS_INSTANTIATE_EXPORT(EigenvalueManager<_1D>)
SOLVERS_INSTANTIATE_EXPORT(EigenvalueManager<_2D>)
SOLVERS_INSTANTIATE_EXPORT(EigenvalueManager<_3D>)

} // end namespace detran


