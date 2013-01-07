//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Matrix.cc
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  Matrix class definition.
 */
//---------------------------------------------------------------------------//

#include "Matrix.hh"
#include "utils/Typedefs.hh"
#include <iostream>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

Matrix::Matrix()
  : d_allocated(false)
{
  /* ... */
}

Matrix::Matrix(const int m, const int n)
  : MatrixBase(m, n)
  , d_allocated(false)
{
  /* ... */
}

Matrix::Matrix(const int m, const int n, const int nnzrow)
  : MatrixBase(m, n)
  , d_allocated(false)
{
  preallocate(nnzrow);
}

Matrix::Matrix(Matrix &A)
  : MatrixBase(A.number_rows(), A.number_columns())
{
  d_nnz = A.number_nonzeros();
  d_rows      = new int[d_m + 1];
  d_columns   = new int[d_nnz];
  d_values    = new double[d_nnz];
  d_diagonals = new int[d_m];

  int size_I = sizeof(d_rows[0]);
  int size_T = sizeof(d_values[0]);
  std::memcpy(d_rows,      A.rows(),      size_I * (d_m + 1) );
  std::memcpy(d_diagonals, A.diagonals(), size_I * d_m       );
  std::memcpy(d_columns,   A.columns(),   size_I * d_nnz     );
  std::memcpy(d_values,    A.values(),    size_T * d_nnz     );

#ifdef CALLOW_ENABLE_PETSC
  PetscErrorCode ierr;
  ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, d_m, d_n,
                                   d_rows, d_columns,
                                   d_values, &d_petsc_matrix);
  Assert(!ierr);
#endif

  d_allocated = true;
  d_is_ready = true;
}

Matrix::~Matrix()
{
  if (d_is_ready)
  {
    delete [] d_values;
    delete [] d_columns;
    delete [] d_rows;
    delete [] d_diagonals;
  }
}

Matrix::SP_matrix
Matrix::Create(const int m, const int n)
{
  SP_matrix p(new Matrix(m, n));
  return p;
}

Matrix::SP_matrix
Matrix::Create(const int m, const int n, const int nnz)
{
  SP_matrix p(new Matrix(m, n, nnz));
  return p;
}

} // end namespace callow

