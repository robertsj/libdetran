//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSolverCreator.cc
 * \brief  EigenSolverCreator 
 * \author Jeremy Roberts
 * \date   Oct 26, 2012
 */
//---------------------------------------------------------------------------//

#include "EigenSolverCreator.hh"
#include "callow/callow_config.hh"
// solvers
#include "PowerIteration.hh"
#include "SlepcSolver.hh"
//
#include <string>

namespace callow
{

/**
  *  Create a solver
  *
  *  @param  tol     tolerance
  *  @param  maxit   maximum number of iterations
  */
EigenSolverCreator::SP_solver
EigenSolverCreator::Create(SP_db db)
{
  SP_solver solver;

  // Set defaults
  std::string solver_type = "power";
  double tol = 1e-5;
  int maxit  = 100;
  int monitor_level = 0;

  // Check database for parameters--easily add new ones here.
  if (db)
  {
    if (db->check("eigen_solver_type"))
      solver_type = db->get<std::string>("eigen_solver_type");
    if (db->check("eigen_solver_tol"))
      tol = db->get<double>("eigen_solver_tol");
    if (db->check("eigen_solver_maxit"))
      maxit = db->get<int>("eigen_solver_maxit");
    if (db->check("eigen_solver_monitor_level"))
      monitor_level = db->get<int>("eigen_solver_monitor_level");
  }

  if (solver_type == "power")
    solver = new PowerIteration(tol, maxit);
  else if (solver_type == "slepc")
  {
#ifdef CALLOW_ENABLE_SLEPC
    solver = new SlepcSolver(tol, maxit);
#else
    THROW("SLEPc solvers not available with this build");
#endif
  }
  else
    THROW("Unsupported solver type: " + solver_type);

  solver->set_monitor_level(monitor_level);
  return solver;
}


} // end namespace callow

//---------------------------------------------------------------------------//
//              end of file EigenSolverCreator.cc
//---------------------------------------------------------------------------//
