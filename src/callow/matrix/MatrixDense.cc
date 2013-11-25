//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixDense.cc
 *  @brief MatrixDense member definitins
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MatrixDense.hh"

namespace callow
{

//----------------------------------------------------------------------------//
MatrixDense::MatrixDense(const int m, const int n, const double v)
  : MatrixBase(m, n)
{
  Require(d_m > 0);
  Require(d_n > 0);

  int N = m * n;
  d_values = new double[N];
  for (int i = 0; i < N; ++i) d_values[i] = v;

#ifdef CALLOW_ENABLE_PETSC
  PetscErrorCode ierr;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, d_n, d_m,
                           d_values, &d_petsc_matrix);
  Assert(!ierr);
#endif

  d_is_ready = true;
}

//----------------------------------------------------------------------------//
MatrixDense::MatrixDense(const MatrixDense &A)
  : MatrixBase(A.number_rows(), A.number_columns())
  , d_values(NULL)
{
  Require(d_m > 0);
  Require(d_n > 0);

  d_values = new double[d_m * d_n];
//  for (int i = 0; i < d_m; ++i)
//    for (int j = 0; j < d_n; ++j)
//      d_values[MDIDX(i, j)] = A(i, j);
  for (int i = 0; i < d_m*d_n; ++i)
    d_values[i] = A[i];
  d_is_ready = true;

#ifdef CALLOW_ENABLE_PETSC
  PetscErrorCode ierr;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, d_n, d_m,
                           d_values, &d_petsc_matrix);
  Assert(!ierr);
#endif
}

//----------------------------------------------------------------------------//
MatrixDense::~MatrixDense()
{
  if (d_is_ready)
  {
//#ifdef CALLOW_ENABLE_PETSC
//  PetscErrorCode ierr;
//  ierr = MatDestroy(&d_petsc_matrix);
//  Assert(!ierr);
//#endif
    delete [] d_values;
  }
}

//----------------------------------------------------------------------------//
MatrixDense::SP_matrix
MatrixDense::Create(const int m, const int n, const double v)
{
  SP_matrix p(new MatrixDense(m, n, v));
  return p;
}


//----------------------------------------------------------------------------//
void MatrixDense::display(bool forceprint) const
{
  Require(d_is_ready);
  printf(" Dense matrix \n");
  printf(" ---------------------------\n");
  printf("      number rows = %5i \n",   d_m);
  printf("   number columns = %5i \n",   d_n);
  printf("\n");
  if ((d_m > 20 || d_n > 20) && !forceprint)
  {
    printf("  *** matrix not printed for m or n > 20 *** \n");
    return;
  }
  for (int i = 0; i < d_m; i++)
  {
    printf(" row  %3i | ", i);
    for (int j = 0; j < d_n; j++)
    {
      double v = d_values[MDIDX(i, j)];
      printf(" %3i (%13.6e)", j, v);
    }
    printf("\n");
  }
  printf("\n");
}

//----------------------------------------------------------------------------//
void MatrixDense::print_matlab(std::string filename) const
{
  Require(d_is_ready);
  FILE * f;
  f = fopen (filename.c_str(), "w");
  for (int i = 0; i < d_m; ++i)
  {
    for (int j = 0; j < d_n; ++j)
    {
      double v = d_values[MDIDX(i, j)];
      fprintf(f, "%8i   %8i    %23.16e \n", i+1, j+1, v);
    }
  }
  fprintf(f, "\n");
  fclose (f);
}

//----------------------------------------------------------------------------//
void MatrixDense::clear()
{
  for (int i = 0; i < d_m; ++i)
  {
    for (int j = 0; j < d_n; ++j)
    {
      d_values[MDIDX(i, j)] = 0.0;
    }
  }
}


} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file MatrixDense.cc
//----------------------------------------------------------------------------//
