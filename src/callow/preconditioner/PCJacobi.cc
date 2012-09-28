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

PCJacobi::PCJacobi(SP_matrix A)
{
  // preconditions
  Require(A);
  Require(A->number_rows() == A->number_columns());

  // create the diagonal vector
  int n = A->number_rows();
  d_P = new Vector(n, 0.0);

  // load values from A and avoid divide by zero
  for (int i = 0; i < n; ++i)
  {
    int d = A->diagonal(i);
    double aii = (*A)[d];
    if (aii == 0.0)
      (*d_P)[i] = 1.0; // or n?
    else
      (*d_P)[i] = 1.0 / aii;  // / aii;
  }
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file PCJacobi.cc
//---------------------------------------------------------------------------//
