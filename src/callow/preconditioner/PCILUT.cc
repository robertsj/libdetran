//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCILUT.cc
 *  @brief PCILUT member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "PCILUT.hh"
#include "SparseRow.hh"

#include <cmath>
#include <algorithm>

namespace callow
{


//----------------------------------------------------------------------------//
PCILUT::PCILUT(SP_matrix A, const size_t p, const double t)
  : Base("PCILUT")
  , d_p(p)
  , d_t(t)
{
  Require(A);
  Require(A->number_rows() == A->number_columns());
  Insist(dynamic_cast<Matrix*>(A.bp()),
    "Need an explicit matrix for use with PCILUT");
  SP_matrixfull B(A);

  d_size = A->number_columns();

  // copy A
  d_P = new Matrix(*B);

  // get the csr structure
  int n = d_P->number_rows();
  double *luval = d_P->values();

  /*
    *    given A, p, and t
    *    initialize L = eye(size(A)) and U = A
    *    row[0:n] = 0
    *    for i = 1:n
    *      row[:] = sparsecopy(a[i, :])
    *      for k = 0:i-1
    *        continue if row[k] == 0
    *        row[k] /= U[k, k]
    *        row[k] = drop_t(row[k], t)
    *        continue if row[k] == 0
    *        for j = k+1:n
    *           row[j] = sparseadd(row[j], -row[k]*U[k, j])
    *        end j
    *      end k
    *      row[0:n] = drop_p(row[0:n], p, t)
    *      L[i, 0:i-1] = sparsecopy(row[0:i-1])
    *      U[i, i:n] = sparsecopy(row[i:n])
    *      row[0:n] = 0
    *    end i
    */

  bool drop_t(SparseRow  &r, int idx, double t)
  {
    if (std::abs(r[idx]) < 0)
    {
      r.delete(idx);
      return false;
    }
    return true;
  }

  bool drop_p(SparseRow  &r, int d, int p, double t)
  {
    // first drop small values
    for (int idx = 0; idx < r.number_nonzeros(); ++idx)
    {
      drop_t(r, idx, t);
    }
    // then keep the p largest on either side of d
    int cL = 0;
    int idx_d = 0;
    for (int idx = 0; idx < r.number_nonzeros(); ++idx)
    {
      int j = r.indices()[idx];
      if (j == d)
        idx_d = idx;
      else if (j < d)
        cL++;
    }
    int cU = r.number_nonzeros() - 1 - cL;
    Ensure(cU >= 0);
    // Now need to find p largest among cU and p largest among
    auto indices = r.indices();
    std::sort(indices.begin())

    return true;
  }

  // P = L + U
  bool keep;
  for (int i = 0; i < n; ++i)
  {
    SparseRow row(*d_P, i);
    for (int k = 0; k < i; ++k)
    {
      int p = row.find(k);
      if (p < 0)
        continue;

      double dv = ((*d_P)[d_P->diagonal(i)]);
      if (dv == 0)
      {
        THROW("ZERO PIVOT IN ILUT");
      }
      row[p] /= dv;
      keep = drop_t(row, p, d_t);
      if (!keep)
        continue;
      SparseRow row_U(*d_P, k);
      row.insert(k+1, n, row_U, -row[k], SparseRow::ADD);
    } // end k
    keep = drop_p(row, i, d_p, d_t);

  }

  delete [] iw;

  // size the working vector
  d_y.resize(d_P->number_rows(), 0.0);

}

//----------------------------------------------------------------------------//
PCILUT::SP_preconditioner
PCILUT::Create(SP_matrix A, const size_t p, const double t)
{
  SP_preconditioner pc(new PCILUT(A, p, t));
  return pc;
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file PCILU0.cc
//----------------------------------------------------------------------------//
