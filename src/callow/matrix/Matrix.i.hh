//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Matrix.i.hh
 *  @brief Matrix inline member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_MATRIX_I_HH_
#define callow_MATRIX_I_HH_

#include "utilities/DBC.hh"
#include <algorithm>
#include <cstring>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>

namespace callow
{

//----------------------------------------------------------------------------//
// PREALLOCATION
//----------------------------------------------------------------------------//

inline void Matrix::preallocate(const int nnzrow)
{
  Require(d_sizes_set);
  Require(!d_allocated);
  Require(nnzrow > 0);
  // set number of nonzeros, preallocate triplets, and initialize counter
  d_aij.resize(d_m, std::vector<triplet_T>(nnzrow));
  d_counter.resize(d_m, 0);
  d_allocated = true;
}


inline void Matrix::preallocate(int* nnzrows)
{
  Require(d_sizes_set);
  Require(!d_allocated);
  // set number of nonzeros, preallocate triplets, and initialize counter
  d_aij.resize(d_m);
  for (int i = 0; i < d_m; i++)
  {
    Require(nnzrows[i] >= 0);
    d_aij[i].resize(nnzrows[i]);
  }
  d_counter.resize(d_m, 0);
  d_allocated = true;
}

//----------------------------------------------------------------------------//
// ASSEMBLING
//----------------------------------------------------------------------------//

/*
 *  We construct in COO format, ie with (i,j,v) triplets.  This makes
 *  adding entries much easier.  Now, we need to order everything
 *  so that the resulting CSR storage has for all i-->j, with j in
 *  order.  The diagonal pointer is also set.  If (i,i,v) doesn't exist,
 *  we insert it, since that's needed for Gauss-Seidel, preconditioning,
 *  and other things.
 */

inline void Matrix::assemble()
{
  Insist(!d_is_ready, "This matrix must not have been assembled already.");
  Insist(d_allocated, "This matrix must be allocated before assembling.");

  typedef std::vector<triplet_T>::iterator it_T;

  // sort the coo structure
  d_nnz = 0;
  for (int i = 0; i < d_m; i++)
  {
//    std::cout << " before sort" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // do the sort
    std::sort(d_aij[i].begin(), d_aij[i].end(), compare_triplet);

//    std::cout << " after sort" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // remove empty entries (if not entered, i,j==-1)
    while ( (d_aij[i].end()-1)->j == -1 ) d_aij[i].pop_back();

//    std::cout << " after pop" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // find and/or insert the diagonal
    Assert(d_aij[i].size());
    if (i < d_n)
    {
      int d = 0;
      while(d < d_aij[i].size())
      {
        if (d_aij[i][d].j == d_aij[i][d].i)
        {
          d = 0;
          break;
        }
        else if (d_aij[i][d].j > d_aij[i][d].i)
        {
          d = d + 1;
          break;
        }
        ++d;
      }
      --d;
      if (d >= 0 && i < d_n)
        d_aij[i].insert(d_aij[i].begin()+d, triplet_T(i, i, 0.0));
    }

//    std::cout << " after diag" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // update the total nonzero count
    d_nnz += d_aij[i].size();
  }

  // allocate
  d_rows = new int[d_m + 1];
  d_columns = new int[d_nnz];
  d_values = new double[d_nnz];
  d_diagonals = new int[d_m];

  // fill the csr storage
  d_rows[0] = 0;
  // pointer to value index
  int p = 0;
  // for all rows
  for (int i = 0; i < d_m; i++)
  {
    // for all columns in the row
    for (size_t j = 0; j < d_aij[i].size(); ++j, ++p)
    {
      d_columns[p] = d_aij[i][j].j;
      d_values[p] = d_aij[i][j].v;
      // store diagonal index
      if (d_columns[p] == i) d_diagonals[i] = p;
    }
    d_rows[i + 1] = p;
    // delete this row of coo storage
    d_aij[i].clear();
  }
  // delete the coo storage and specify the matrix is set to use
  d_aij.clear();

#ifdef CALLOW_ENABLE_PETSC
  PetscErrorCode ierr;
  ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, d_m, d_n,
                                   d_rows, d_columns,
                                   d_values, &d_petsc_matrix);
  Assert(!ierr);
#endif
  d_is_ready = true;
}

//----------------------------------------------------------------------------//
// INDEXING
//----------------------------------------------------------------------------//

inline int Matrix::start(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  return d_rows[i];
}

inline int Matrix::diagonal(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  return d_diagonals[i];
}


inline int Matrix::end(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  return d_rows[i + 1];
}


inline int Matrix::column(const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 && p < d_nnz);
  return d_columns[p];
}


inline double Matrix::operator[](const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 && p < d_nnz);
  return d_values[p];
}

//----------------------------------------------------------------------------//
// ACCESS
//----------------------------------------------------------------------------//

inline double Matrix::operator()(const int i, const int j) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  Require(j >= 0 && j < d_n);
  // loop through elements in this row
  if (i == j)
    return d_values[d_diagonals[i]];
  else if (i < j)
  {
    for (int k = d_rows[i]; k < d_diagonals[i]; k++)
      if (d_columns[k] == j) return d_values[k];
  }
  else
  {
    for (int k = d_diagonals[i] + 1; k < d_rows[i + 1]; k++)
      if (d_columns[k] == j) return d_values[k];
  }
  // otherwise, this is a zeros
  return 0.0;
}

//----------------------------------------------------------------------------//
// MULTIPLY
//----------------------------------------------------------------------------//

inline void Matrix::multiply(const Vector &x, Vector &y)
{
  Require(d_is_ready);
  Require(x.size() == d_n);
  Require(y.size() == d_m);
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMult(d_petsc_matrix, const_cast<Vector* >(&x)->petsc_vector(), y.petsc_vector());
#else
  // clear the output vector
  y.scale(0);
  // row pointer
  int p = 0;
  // for all rows
  #pragma omp for private(p)
  for (int i = 0; i < d_m; i++)
  {
    double temp = y[i];
    // for all columns
    for (p = d_rows[i]; p < d_rows[i + 1]; p++)
    {
      int j = d_columns[p];
      temp += x[j] * d_values[p];
    }
    y[i] = temp;
  }
#endif
}

// \todo good threading option needed

inline void Matrix::multiply_transpose(const Vector &x, Vector &y)
{
  Require(d_is_ready);
  Require(x.size() == d_m);
  Require(y.size() == d_n);
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMultTranspose(d_petsc_matrix,
                   const_cast<Vector* >(&x)->petsc_vector(),
                   y.petsc_vector());
#else
  // clear the output vector
  y.scale(0);
  // for all rows (now columns)
  for (int i = 0; i < d_m; i++)
  {
    // for all columns (now rows)
    for (int p = d_rows[i]; p < d_rows[i + 1]; p++)
    {
      int j = d_columns[p];
      y[j] += x[i] * d_values[p];
    }
  }
#endif
}

//----------------------------------------------------------------------------//
// INSERTING VALUES
//----------------------------------------------------------------------------//


inline bool Matrix::insert(int i, int j, double v, const int type)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 && i < d_m);
  Require(j >= 0 && j < d_n);
  // add to a current entry if found
  if (type == ADD)
  {
    for (size_t p = 0; p < d_aij[i].size(); ++p)
    {
      if (d_aij[i][p].j == j)
      {
        d_aij[i][p].v += v;
        return true;
      }
    }
  }
  // otherwise it didn't exist and return if storage unavailable
  if (d_counter[i] >= d_aij[i].size()) return false;
  // otherwise, add the entry
  d_aij[i][d_counter[i]].i = i;
  d_aij[i][d_counter[i]].j = j;
  d_aij[i][d_counter[i]].v = v;
  ++d_counter[i];
  return true;
}


inline bool Matrix::insert(int i, int *j, double *v, int n, const int type)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 && i < d_m);
  if (type == ADD)
  {
    // loop through all columns given
    for (int jj = 0; jj < n; ++jj)
    {
      bool flag = false;
      for (size_t p = 0; p < d_aij[i].size(); ++p)
      {
        if (d_aij[i][p].j == j[jj])
        {
          d_aij[i][p].v += v[jj];
          flag = true;
          break;
        }
      }
      if (!flag) // insert if not present
        flag = insert(i, j[jj], v[jj], INSERT);
      if (!flag) return false; // fail if we can't insert
    }
    return true;
  }
  // otherwise it didn't exist and return if storage unavailable
  // return if storage unavailable
  if (d_counter[i] + n >  d_aij[i].size()) return false;
  // otherwise, add the entries
  for (int jj = 0; jj < n; ++jj)
  {
    Require(j[jj] >= 0 && j[jj] < d_n);
    d_aij[i][d_counter[i]].i = i;
    d_aij[i][d_counter[i]].j = j[jj];
    d_aij[i][d_counter[i]].v = v[jj];
    ++d_counter[i];
  }
  return true;
}


inline bool Matrix::insert(int *i, int j, double *v, int n, const int type)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(j >= 0 && j < d_n);
  Insist(type == INSERT, "Cannot ADD by column");
  // return if storage unavailable
  for (int ii = 0; ii < n; ++ii)
  {
    Require(i[ii] >= 0 && i[ii] < d_m);
    if (d_counter[i[ii]] + 1 > d_aij[i[ii]].size())
      return false;
  }
  // otherwise, add the entries
  for (int ii = 0; ii < n; ++ii)
  {
    d_aij[i[ii]][d_counter[i[ii]]].i = i[ii];
    d_aij[i[ii]][d_counter[i[ii]]].j = j;
    d_aij[i[ii]][d_counter[i[ii]]].v = v[ii];
    ++d_counter[i[ii]];
  }
  return true;
}


inline bool Matrix::insert(int *i, int *j, double *v, int n, const int type)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Insist(type == INSERT, "Cannot ADD by row, column, value triplet")
  // return if storage unavailable
  // \todo this assumes one entry per row---fix
  for (int k = 0; k < n; ++k)
  {
    Require(i[k] >= 0 && i[k] < d_m);
    Require(j[k] >= 0 && j[k] < d_n);
    if (d_counter[i[k]] + 1 > d_aij[i[k]].size())
      return false;
  }
  // otherwise, add the entries
  for (int k = 0; k < n; ++k)
  {
    d_aij[i[k]][d_counter[i[k]]].i = i[k];
    d_aij[i[k]][d_counter[i[k]]].j = j[k];
    d_aij[i[k]][d_counter[i[k]]].v = v[k];
    ++d_counter[i[k]];
  }
  return true;
}

} // end namespace callow

#endif /* callow_MATRIX_I_HH_ */

//----------------------------------------------------------------------------//
//              end of Matrix.i.hh
//----------------------------------------------------------------------------//
