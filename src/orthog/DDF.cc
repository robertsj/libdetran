//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DDF.cc
 *  @brief DDF member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "DDF.hh"
#include "utilities/DBC.hh"

namespace detran_orthog
{

//----------------------------------------------------------------------------//
DDF::DDF(const Parameters &p)
: OrthogonalBasis(p)
{
 Insist(d_order + 1 == d_size,
   "It makes no sense to use DDF for incomplete expansions!");

 // Allocate the basis matrix and fill diagonal.
 d_basis = new callow::MatrixDense(d_size, d_size, 0.0);
 for (size_t i = 0; i < d_order + 1; ++i)
   (*d_basis)(i, i) = 1.0;
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file DDF.cc
//----------------------------------------------------------------------------//



