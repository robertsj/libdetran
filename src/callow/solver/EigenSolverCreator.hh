//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  EigenSolverCreator.hh
 *  @brief EigenSolverCreator class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef callow_EIGENSOLVERCREATOR_HH_
#define callow_EIGENSOLVERCREATOR_HH_

#include "EigenSolver.hh"

namespace callow
{

/**
 *  @class EigenSolverCreator
 *  @brief Creates an eigensolver
 */
class CALLOW_EXPORT EigenSolverCreator
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

//----------------------------------------------------------------------------//
//              end of file EigenSolverCreator.hh
//----------------------------------------------------------------------------//
