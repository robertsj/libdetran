//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearSolverCreator.hh
 *  @brief  LinearSolverCreator
 *  @author Jeremy Roberts
 *  @date   Sep 21, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_SOLVERCREATER_HH_
#define callow_SOLVERCREATER_HH_

#include "LinearSolver.hh"
#include "utilities/InputDB.hh"

#include <string>

namespace callow
{

/**
 *  @class LinearSolverCreator
 *  @brief Creates a
 */
class CALLOW_EXPORT LinearSolverCreator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver::SP_solver               SP_solver;
  typedef detran_utilities::InputDB::SP_input   SP_db;

  //-------------------------------------------------------------------------//
  // PUBLIC METHODS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Create a linear solver
   *  @param  db  Pointer to parameter database
   */
  static SP_solver Create(SP_db db = SP_db(0));

};

} // end namespace callow

#endif // callow_SOLVERCREATOR_HH_

//---------------------------------------------------------------------------//
//              end of file SolverCreator.hh
//---------------------------------------------------------------------------//
