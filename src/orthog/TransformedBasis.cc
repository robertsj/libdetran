//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TransformedBasis.cc
 *  @brief TransformedBasis member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#include "TransformedBasis.hh"
#include "DCT.hh"
#include "utilities/DBC.hh"

namespace detran_orthog
{

using namespace callow;

//----------------------------------------------------------------------------//
TransformedBasis::TransformedBasis(const Parameters &p)
  : OrthogonalBasis(p)
{
  Require(p.x.size() == d_size);

  d_orthonormal = true;

  Vector x(p.x);
  double norm_x = x.norm(callow::L2);
  Assert(norm_x != 0.0);

  // Allocate the base DCT matrix
  DCT dct(p);
  MatrixDense &P = *dct.basis();
  d_basis = new callow::MatrixDense(P);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // Set the zeroth order
  for (int i = 0; i < d_size; ++i)
  {
    (*d_basis)(i, 0) *= x[i];
  }
  Vector B_0(d_size, &(*d_basis)(0, 0));
  B_0.scale(1.0/B_0.norm(L2));

  // Orthogonalize via Gram-Schmidt
  for (int l = 1; l <= d_order; ++l)
  {
    Vector B_l(d_size, &(*d_basis)(0, l));
    for (int m = 0; m < l; ++m)
    {
      Vector B_m(d_size, &(*d_basis)(0, m));
      double c = -B_m.dot(B_l);
      B_l.add_a_times_x(c, B_m);
    }
    B_l.scale(1.0/B_l.norm(L2));
  }

  compute_a();
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file TransformedBasis.cc
//----------------------------------------------------------------------------//




