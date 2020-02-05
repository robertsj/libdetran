//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  UserBasis.cc
 *  @brief UserBasis member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "UserBasis.hh"
#include "utilities/DBC.hh"

namespace detran_orthog
{

//----------------------------------------------------------------------------//
UserBasis::UserBasis(const Parameters &p)
  : OrthogonalBasis(p)
{
  d_orthonormal = true;

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_size, d_order + 1, 0.0);


  // Build the basis
  for (size_t l = 0; l <= d_order; ++l)
  {
    vec_dbl USR_funct = p.db->get<vec_dbl>("vec"+AsString(l));
    for (size_t i = 0; i < d_size; ++i)
    {
     (*d_basis)(i, l) = USR_funct[i];
    }
  }
  compute_a();
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file UserBasis.cc
//----------------------------------------------------------------------------//
