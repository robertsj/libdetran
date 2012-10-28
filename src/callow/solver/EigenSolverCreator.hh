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

namespace callow
{

/**
 *  \class EigenSolverCreator
 *  \brief Creates an eigensolver
 */
class EigenSolverCreator
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<EigenSolver>   SP_solver;
  typedef detran_utilities::InputDB::SP_input SP_db;

  //-------------------------------------------------------------------------//
  // PUBLIC METHODS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Create an eigensolver
   *  @param  db  Pointer to parameter database
   */
  static SP_solver Create(SP_db db = SP_db(0));

};

} // end namespace callow


#endif /* callow_EIGENSOLVERCREATOR_HH_ */
