//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ReflectiveSolver.cc
 * \brief  ReflectiveSolver 
 * \author Jeremy Roberts
 * \date   Nov 2, 2012
 */
//---------------------------------------------------------------------------//

#include "ReflectiveSolver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
ReflectiveSolver<D>::ReflectiveSolver(SP_state        state,
                                      SP_boundary     boundary,
                                      SP_sweeper      sweeper,
                                      SP_sweepsource  source)
  : d_state(state)
  , d_boundary(boundary)
  , d_sweeper(sweeper)
  , d_source(source)
{
  // Preconditions
  Require(d_state);
  Require(d_boundary);
  Require(d_sweeper);
  Require(d_source);

  // Get input
  d_input = d_state->get_input();
  Require(d_input);

}

//


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file ReflectiveSolver.cc
//---------------------------------------------------------------------------//
