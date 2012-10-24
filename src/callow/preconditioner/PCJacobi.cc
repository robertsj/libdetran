//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCJacobi.cc
 * \brief  PCJacobi 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#include "PCJacobi.hh"

namespace callow
{

//---------------------------------------------------------------------------//
PCJacobi::PCJacobi(SP_matrix A)
{
  // preconditions
  Require(A);
  Require(A->number_rows() == A->number_columns());
  Insist(dynamic_cast<Matrix*>(A.bp()),
    "Need an explicit matrix for use with PCILU0");
  SP_matrixfull B(A);

  // create the diagonal vector
  int n = B->number_rows();
  d_P = new Vector(n, 0.0);

  // load values from A and avoid divide by zero
  for (int i = 0; i < n; ++i)
  {
    int d = B->diagonal(i);
    double aii = (*B)[d];
    if (aii == 0.0)
      (*d_P)[i] = 1.0; // or n?
    else
      (*d_P)[i] = 1.0 / aii;  // / aii;
  }
}

//---------------------------------------------------------------------------//
PCJacobi::SP_preconditioner PCJacobi::Create(SP_matrix A)
{
  SP_preconditioner p(new PCJacobi(A));
  return p;
}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file PCJacobi.cc
//---------------------------------------------------------------------------//
