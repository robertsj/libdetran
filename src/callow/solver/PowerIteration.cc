//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIteration.cc
 * \brief  PowerIteration 
 * \author Jeremy Roberts
 * \date   Sep 24, 2012
 */
//---------------------------------------------------------------------------//

#include "PowerIteration.hh"

namespace callow
{

PowerIteration::PowerIteration(const double    tol,
                                  const int       maxit)
  : Base(tol, maxit, "power")
{
  /* ... */
}

} // end namespace callow

//---------------------------------------------------------------------------//
//              end of file PowerIteration.cc
//---------------------------------------------------------------------------//
