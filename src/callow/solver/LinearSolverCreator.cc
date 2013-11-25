//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearSolverCreator.cc
 *  @brief  LinearSolverCreator
 *  @author Jeremy Roberts
 *  @date   Sep 21, 2012
 */
//---------------------------------------------------------------------------//

#include "LinearSolverCreator.hh"
// solvers
#include "Richardson.hh"
#include "Jacobi.hh"
#include "GaussSeidel.hh"
#include "GMRES.hh"
#include "PetscSolver.hh"

namespace callow
{

//---------------------------------------------------------------------------//
LinearSolverCreator::SP_solver
LinearSolverCreator::Create(SP_db db)
{
  SP_solver solver;

  // Set defaults
  std::string solver_type = "gmres";
  double atol = 1e-8;
  double rtol = 1e-8;
  int maxit   = 100;
  int monitor_level = 0;
  bool monitor_diverge = false;
  double omega = 1.0;
  int restart = 30;

  if (db)
  {
    if (db->check("linear_solver_type"))
      solver_type = db->get<std::string>("linear_solver_type");
    if (db->check("linear_solver_atol"))
      atol = db->get<double>("linear_solver_atol");
    if (db->check("linear_solver_rtol"))
      rtol = db->get<double>("linear_solver_rtol");
    if (db->check("linear_solver_maxit"))
      maxit = db->get<int>("linear_solver_maxit");
    if (db->check("linear_solver_monitor_level"))
      monitor_level = db->get<int>("linear_solver_monitor_level");
    if (db->check("linear_solver_monitor_diverge"))
      monitor_diverge = db->get<int>("linear_solver_monitor_diverge");
    if (solver_type == "richardson" &&
        db->check("linear_solver_richardson_omega"))
    {
      omega = db->get<double>("linear_solver_richardson_omega");
    }
    if (solver_type == "sor" &&
        db->check("linear_solver_sor_omega"))
    {
      omega = db->get<double>("linear_solver_sor_omega");
    }
    if (solver_type == "gmres" &&
        db->check("linear_solver_gmres_restart"))
    {
      restart = db->get<int>("linear_solver_gmres_restart");
    }
  }

//  std::cout << " CALLOW:" << std::endl;
//  std::cout << "  solver_type = " << solver_type << std::endl;
//  std::cout << "         atol = " << atol << std::endl;
//  std::cout << "         rtol = " << rtol << std::endl;
//  std::cout << "        maxit = " << maxit << std::endl;

  //---------------------------------------------------------------------------//
  if (solver_type == "richardson")
  {
    solver = new Richardson(atol, rtol, maxit, omega);
  }

  //---------------------------------------------------------------------------//
  else if (solver_type == "jacobi")
  {
    solver = new Jacobi(atol, rtol, maxit);
  }

  //---------------------------------------------------------------------------//
  else if (solver_type == "gauss-seidel")
  {
    solver = new GaussSeidel(atol, rtol, maxit);
  }

  //---------------------------------------------------------------------------//
  else if (solver_type == "sor")
  {
    solver = new GaussSeidel(atol, rtol, maxit, omega);
  }

  //---------------------------------------------------------------------------//
  else if (solver_type == "gmres")
  {
    solver = new GMRES(atol, rtol, maxit, restart);
  }

  //---------------------------------------------------------------------------//
  else if (solver_type == "petsc")
  {
#ifdef CALLOW_ENABLE_PETSC
    solver = new PetscSolver(atol, rtol, maxit);
#else
    THROW("PETSc solvers not available with this build");
#endif
  }

  //---------------------------------------------------------------------------//
  else
  {
    THROW("Unsupported solver type: " + solver_type);
  }

  solver->set_monitor_diverge(monitor_diverge);
  solver->set_monitor_level(monitor_level);

  return solver;
}

} // end namespace callow

//---------------------------------------------------------------------------//
//              end of file LinearSolverCreator.cc
//---------------------------------------------------------------------------//
