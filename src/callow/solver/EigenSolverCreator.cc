//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenSolverCreator.cc
 *  @brief EigenSolverCreator
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenSolverCreator.hh"
#include "callow/callow_config.hh"
// solvers
#include "PowerIteration.hh"
#include "SlepcSolver.hh"
#include "Davidson.hh"
#include "Eispack.hh"
//
#include <string>

namespace callow
{

//----------------------------------------------------------------------------//
EigenSolverCreator::SP_solver
EigenSolverCreator::Create(SP_db db)
{
  SP_solver solver;

  // Set defaults
  std::string solver_type = "power";
  double tol = 1e-5;
  int maxit = 100;
  int monitor_level = 0;
  int number_values = 1;
  int subspace_size = 20;

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
    if (db->check("eigen_number_values"))
      number_values = db->get<int>("eigen_number_values");
    // note, this is ignored unless we use slepc
  }

  if (solver_type == "power")
    solver = new PowerIteration(tol, maxit);
  else if (solver_type == "slepc")
  {
#ifdef CALLOW_ENABLE_SLEPC
    std::string slepc_type = "krylovschur";
    if (db->check("eigen_solver_slepc_type"))
      slepc_type = db->get<std::string>("eigen_solver_slepc_type");
    solver = new SlepcSolver(slepc_type, tol, maxit, number_values);
#else
    THROW("SLEPc solvers not available with this build");
#endif
  }
  else if (solver_type == "gd")
  {
    if (db->check("eigen_solver_subspace_size"))
      subspace_size = db->get<int>("eigen_solver_subspace_size");
    solver = new Davidson(tol, maxit, subspace_size);
  }
  else if (solver_type == "eispack")
  {
    if (db->check("eigen_solver_subspace_size"))
      subspace_size = db->get<int>("eigen_solver_subspace_size");
    solver = new Eispack(tol, maxit);
  }
  else
    THROW("Unsupported solver type: " + solver_type);

  solver->set_monitor_level(monitor_level);
  return solver;
}


} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file EigenSolverCreator.cc
//----------------------------------------------------------------------------//
