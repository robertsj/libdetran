//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PowerIteration.cc
 *  @brief PowerIteration member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "PowerIteration.hh"

namespace callow
{

//----------------------------------------------------------------------------//
PowerIteration::PowerIteration(const double tol,
                               const int    maxit)
  : Base(tol, maxit, "power")
{
  /* ... */
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file PowerIteration.cc
//----------------------------------------------------------------------------//
