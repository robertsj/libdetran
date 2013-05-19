//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  JFNK.hh
 *  @brief JFNK class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_JFNK_HH_
#define detran_JFNK_HH_

#include "Eigensolver.hh"
#include "EnergyIndependentEigenOperator.hh"

namespace detran
{

/**
 *  @class JFNK
 *  @brief Solves the eigenvalue problem via Jacobian-Free Newton-Krylov
 *
 *  The eigenvalue problem can be cast in the nonlinear form
 *
 *    f = | (A-k*I)*fd    | = 0
 *        | 1/2 - fd'*fd |
 *
 *  where A is the monoenergetic eigenoperator.
 *
 *  The Jacobian is
 *
 *    J = | (A-k*I)    fd | = 0
 *        |  fd'       0  |
 *
 */

} // end namespace detran

#endif // detran_JFNK_HH_

//----------------------------------------------------------------------------//
//              end of file JFNK.hh
//----------------------------------------------------------------------------//
