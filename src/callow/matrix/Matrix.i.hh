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
// ACCESS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
inline int Matrix::start(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  return d_rows[i];
}

//----------------------------------------------------------------------------//
inline int Matrix::diagonal(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  return d_diagonals[i];
}

//----------------------------------------------------------------------------//
inline int Matrix::end(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  return d_rows[i + 1];
}

//----------------------------------------------------------------------------//
inline int Matrix::column(const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 && p < d_nnz);
  return d_columns[p];
}

//----------------------------------------------------------------------------//
inline double Matrix::operator[](const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 && p < d_nnz);
  return d_values[p];
}

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

//----------------------------------------------------------------------------//
inline void Matrix::multiply(const Vector &x, Vector &y)
{
  Require(d_is_ready);
  Require(x.size() == d_n);
  Require(y.size() == d_m);
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMult(d_petsc_matrix,
          const_cast<Vector* >(&x)->petsc_vector(), y.petsc_vector());
#else
  // clear the output vector
  y.scale(0);
  // row pointer
  int p = 0;
  // for all rows
  #pragma omp for schedule(static) private(p)
  for (int i = 0; i < d_m; ++i)
  {
    double temp = y[i];
    // for all columns
    for (p = d_rows[i]; p < d_rows[i + 1]; ++p)
    {
      int j = d_columns[p];
      temp += x[j] * d_values[p];
    }
    y[i] = temp;
  }
#endif
}

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
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
