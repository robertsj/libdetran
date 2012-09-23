//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCJacobi.i.hh
 * \brief  PCJacobi.i 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PCJACOBI_I_HH_
#define PCJACOBI_I_HH_

namespace callow
{

template <class T>
PCJacobi<T>::PCJacobi(SP_matrix A)
{
  // preconditions
  Require(A);
  Require(A->number_rows() == A->number_columns());

  // create the diagonal vector
  int n = A->number_rows();
  d_P = new Vector<T>(n, 0.0);

  // load values from A and avoid divide by zero
  for (int i = 0; i < n; ++i)
  {
    int d = A->diagonal(i);
    T aii = (*A)[d];
    if (aii == 0.0)
      (*d_P)[i] = 1.0; // or n?
    else
      (*d_P)[i] = 1.0 / aii;  // / aii;
  }
  d_P->print_matlab("pcjacobi.out");
}

template <class T>
void PCJacobi<T>::apply(Vector<T> &b, Vector<T> &x)
{
  // preconditions
  Require(x.size() == d_P->size());

  // apply x = inv(P)*b
  x.copy(b);
  x.multiply((*d_P));
}


} // end namespace detran

#endif // PCJACOBI_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file PCJacobi.i.hh
//---------------------------------------------------------------------------//
