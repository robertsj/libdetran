//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCILU0.i.hh
 * \brief  PCILU0.i 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_PCILU0_I_HH_
#define callow_PCILU0_I_HH_

#include <iostream>

namespace callow
{

template <class T>
PCILU0<T>::PCILU0(SP_matrix A)
{
  // preconditions
  Require(A);
  Require(A->number_rows() == A->number_columns());

  // copy A
  d_P = new Matrix<T>(*A);

  using std::cout;
  using std::endl;

  /* the following mostly follows the algorithm of
   * saad in ch 10, which is basically as follows:
   *
   *  for i = 1, m
   *    for k = 0, i-1
   *      if (A(i,k) > 0) A(i,k) = A(i,k)/A(i,i)
   *      for j = k+1, n
   *        if (A(i,j) > 0) A(i,j) = A(i,j) - A(i,j)*A(k,j)
   *
   */

  // working array
  int* iw = new int[A->number_columns()];
  for (int k = 0; k < A->number_columns(); ++k) iw[k] = -1;


  // get the csr structure
  int n = d_P->number_rows();
  T   *luval = d_P->values();

  // loop over rows
  for (int i = 1; i < n; ++i)
  {

    // pre-store the column pointers for this row.  if
    // the column isn't present, the value remains -1
    for (int p = d_P->start(i); p < d_P->diagonal(i); ++p)
      iw[d_P->column(p)] = p;

    // loop through the columns
    for (int p = d_P->start(i); p < d_P->diagonal(i); ++p)
    {
      // column index of row i
      int k = d_P->column(p);

      // compute row multiplier (aik = aik / akk)
      T val = luval[p] / luval[d_P->diagonal(k)];
      luval[p] = val;
      // for *row* k, loop over the columns j above diagonal
      for (int q = d_P->diagonal(k) + 1; q < d_P->end(k); ++q)
      {
        int j = d_P->column(q);
        int pp = iw[j];
        if (pp != -1) luval[pp] -= val * luval[q];
      }
    }
    int d = d_P->diagonal(i);
    if (luval[d] == 0.0)
    {
      THROW("ZERO PIVOT IN ILU0");
    }

    // reset
    for (int p = d_P->start(i); p < d_P->diagonal(i); ++p)
      iw[d_P->column(p)] = -1;
  }

  delete [] iw;

  // size the working vector
  d_y.resize(d_P->number_rows(), 0.0);

}

template <class T>
void PCILU0<T>::apply(Vector<T> &b, Vector<T> &x)
{

  // solve LUx = x --> x = inv(U)*inv(L)*x

  // forward substitution
  //   for i = 0:m-1
  //     x[i] = 1/L[i,i] * ( b[i] - sum(k=0:i-1, L[i,k]*y[k]) )
  // but note that in our ILU(0) scheme, L is *unit* lower triangle,
  // meaning L has ones on the diagonal (whereas U does not)
  for (int i = 0; i < d_P->number_rows(); ++i)
  {
    // start index
    int s = d_P->start(i);
    // diagonal index
    int d = d_P->diagonal(i);
    // invert row
    d_y[i] = b[i];
    for (int p = s; p < d; ++p)
    {
      // column index
      int c = d_P->column(p);
      d_y[i] -= d_P->values()[p] * d_y[c];
    }
  }

  // backward substitution
  //   for i = m-1:0
  //     y[i] = 1/U[i,i] * ( b[i] - sum(k=i+1:m-1, U[i,k]*y[k]) )
  for (int i = d_P->number_rows() - 1; i >= 0; --i)
  {
    // diagonal index
    int d = d_P->diagonal(i);
    // end index
    int e = d_P->end(i);
    // invert row
    x[i] = d_y[i];
    for (int p = d; p < e; ++p)
    {
      // column index
      int c = d_P->column(p);
      x[i] -= d_P->values()[p] * x[c];
    }
    x[i] /= d_P->values()[d];
  }

}

} // end namespace detran

#endif // callow_PCILU0_I_HH_

//---------------------------------------------------------------------------//
//              end of file PCILU0.i.hh
//---------------------------------------------------------------------------//
