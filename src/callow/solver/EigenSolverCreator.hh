//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EigenSolverCreator.hh
 *  @author robertsj
 *  @date   Sep 25, 2012
 *  @brief  EigenSolverCreator class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_EIGENSOLVERCREATOR_HH_
#define callow_EIGENSOLVERCREATOR_HH_

#include "EigenSolver.hh"
// solvers
#include "PowerIteration.hh"
#include "SlepcSolver.hh"
//
#include <string>

namespace callow
{

/**
 *  \class EigenSolverCreator
 *  \brief Creates an eigensolver
 */
template <class T>
class EigenSolverCreator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<EigenSolver<T> >   SP_solver;

  //-------------------------------------------------------------------------//
  // PUBLIC METHODS
  //-------------------------------------------------------------------------//

  /**
   *  Create a solver
   *
   *  @param  tol     tolerance
   *  @param  maxit   maximum number of iterations
   */
  static SP_solver Create(std::string  type,
                          const double tol = 1e-5,
                          const int    maxit = 100)
  {
    SP_solver solver;

    if (type == "power")
      solver = new PowerIteration<T>(tol, maxit);
    else if (type == "slepc")
    {
#ifdef CALLOW_ENABLE_PETSC
      solver = new SlepcSolver(tol,maxit);
#else
      THROW("SLEPc solvers not available with this build");
#endif
    }
    else
      THROW("Unsupported solver type requested");

    return solver;
  }

};

} // end namespace callow


#endif /* callow_EIGENSOLVERCREATOR_HH_ */
