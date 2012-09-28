//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolverCreator.hh
 * \brief  LinearSolverCreator
 * \author Jeremy Roberts
 * \date   Sep 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_SOLVERCREATER_HH_
#define callow_SOLVERCREATER_HH_

#include "LinearSolver.hh"
// solvers
#include "Richardson.hh"
#include "Jacobi.hh"
#include "GaussSeidel.hh"
#include "GMRES.hh"
#include "PetscSolver.hh"

#include <string>

namespace callow
{

/**
 *  @class LinearSolverCreator
 *  @brief Creates a
 */
class LinearSolverCreator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef typename LinearSolver::SP_solver   SP_solver;

  //-------------------------------------------------------------------------//
  // PUBLIC METHODS
  //-------------------------------------------------------------------------//

  /**
   *  Create a solver
   *
   *  @param  atol    absolute tolerance
   *  @param  rtol    relative tolerance
   *  @param  maxit   maximum number of iterations
   */
  static SP_solver Create(std::string  type,
                          const double atol,
                          const double rtol,
                          const int    maxit = 1000)
  {
    SP_solver solver;

    if (type == "richardson")
      solver = new Richardson(atol, rtol, maxit);
    else if (type == "jacobi")
      solver = new Jacobi(atol, rtol, maxit);
    else if (type == "gauss-seidel")
      solver = new GaussSeidel(atol, rtol, maxit);
    else if (type == "gmres")
      solver = new GMRES(atol, rtol, maxit);
    else if (type == "petsc")
    {
#ifdef CALLOW_ENABLE_PETSC
      solver = new PetscSolver(atol, rtol, maxit);
#else
      THROW("PETSc solvers not available with this build");
#endif
    }
    else
      THROW("Unsupported solver type requested");

    return solver;
  }

};

} // end namespace callow

#endif // callow_SOLVERCREATOR_HH_

//---------------------------------------------------------------------------//
//              end of file SolverCreator.hh
//---------------------------------------------------------------------------//
