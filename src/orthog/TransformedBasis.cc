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

  // Allocate the base matrix
  Insist(p.transformed_key != "trans",
    "Cannot transform a transformed basis.");
  SP_basis original = OrthogonalBasis::Create(p.transformed_key, p);
  const MatrixDense &P = *(original->basis());
  d_basis = new callow::MatrixDense(P);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // Set the zeroth order, located in the abscissa vector for now
  for (int i = 0; i < d_size; ++i)
  {
    (*d_basis)(i, 0) *= x[i];
  }
  Vector B_0(d_size, &(*d_basis)(0, 0));
  B_0.scale(1.0/B_0.norm(L2));

  // Orthogonalize via Gram-Schmidt
  for (int l = 1; l <= d_order; ++l)
  {
    // Copy the next order basis from the original basis
    Vector B_l(d_size, &(*d_basis)(0, l));
    // Optionally multiply this pointwise by the zeroth order
    if (p.transformed_option == 1)
    {
      for (int i = 0; i < d_size; ++i)
        B_l[i] *= B_0[i];
    }
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




