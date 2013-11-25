//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCILU0.cc
 *  @brief PCILU0 member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "PCILU0.hh"

namespace callow
{

//----------------------------------------------------------------------------//
PCILU0::PCILU0(SP_matrix A)
  : Base("PCILU0")
{
  Require(A);
  Require(A->number_rows() == A->number_columns());
  Insist(dynamic_cast<Matrix*>(A.bp()),
    "Need an explicit matrix for use with PCILU0");
  SP_matrixfull B(A);

  d_size = A->number_columns();

  // copy A
  d_P = new Matrix(*B);

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
  double *luval = d_P->values();

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
      double val = luval[p] / luval[d_P->diagonal(k)];
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

//----------------------------------------------------------------------------//
PCILU0::SP_preconditioner PCILU0::Create(SP_matrix A)
{
  SP_preconditioner p(new PCILU0(A));
  return p;
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file PCILU0.cc
//----------------------------------------------------------------------------//
