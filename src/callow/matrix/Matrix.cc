//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Matrix.cc
 *  @brief Matrix member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Matrix.hh"
#include "utils/Typedefs.hh"
#include <iostream>

namespace callow
{

//----------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
Matrix::Matrix()
  : d_values(NULL)
  , d_columns(NULL)
  , d_rows(NULL)
  , d_diagonals(NULL)
  , d_nnz(0)
  , d_allocated(false)
{
  /* ... */
}

//----------------------------------------------------------------------------//
Matrix::Matrix(const int m, const int n)
  : MatrixBase(m, n)
  , d_values(NULL)
  , d_columns(NULL)
  , d_rows(NULL)
  , d_diagonals(NULL)
  , d_nnz(0)
  , d_allocated(false)
{
  /* ... */
}

//----------------------------------------------------------------------------//
Matrix::Matrix(const int m, const int n, const int nnzrow)
  : MatrixBase(m, n)
  , d_allocated(false)
{
  preallocate(nnzrow);
}

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
// IO
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
void Matrix::display(bool forceprint) const
{
  Require(d_is_ready);
  printf(" CSR matrix \n");
  printf(" ---------------------------\n");
  printf("      number rows = %5i \n",   d_m);
  printf("   number columns = %5i \n",   d_n);
  printf("      stored size = %5i \n\n", d_nnz);
  printf("\n");
  if ((d_m > 20 || d_n > 20) && !forceprint)
  {
    printf("  *** matrix not printed for m or n > 20 *** \n");
    return;
  }
  for (int i = 0; i < d_m; i++)
  {
    printf(" row  %3i | ", i);
    for (int p = d_rows[i]; p < d_rows[i + 1]; p++)
    {
      int j = d_columns[p];
      double v = d_values[p];
      printf(" %3i (%13.6e)", j, v);
    }
    printf("\n");
  }
  printf("\n");
}

//----------------------------------------------------------------------------//
inline void Matrix::print_matlab(std::string filename) const
{
  Require(d_is_ready);
  FILE * f;
  f = fopen (filename.c_str(), "w");
  for (int i = 0; i < d_m; i++)
  {
    for (int p = d_rows[i]; p < d_rows[i + 1]; p++)
    {
      int j = d_columns[p];
      double v = d_values[p];
      fprintf(f, "%8i   %8i    %23.16e \n", i+1, j+1, v);
    }
  }
  fprintf(f, "\n");
  fclose (f);
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of Matrix.cc
//----------------------------------------------------------------------------//
