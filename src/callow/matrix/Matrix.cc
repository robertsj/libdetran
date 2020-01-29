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
Matrix::SP_matrix Matrix::Create(const int m, const int n)
{
  SP_matrix p(new Matrix(m, n));
  return p;
}

//----------------------------------------------------------------------------//
Matrix::SP_matrix Matrix::Create(const int m, const int n, const int nnz)
{
  SP_matrix p(new Matrix(m, n, nnz));
  return p;
}

//----------------------------------------------------------------------------//
void Matrix::preallocate(const int nnzrow)
{
  Require(d_sizes_set);
  Require(!d_allocated);
  Require(nnzrow > 0);
  // set number of nonzeros, preallocate triplets, and initialize counter
  d_aij.resize(d_m, std::vector<triplet_T>(nnzrow));
  d_counter.resize(d_m, 0);
  d_allocated = true;
}

//----------------------------------------------------------------------------//
void Matrix::preallocate(int* nnzrows)
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

//#define db_print_row(i, c)                  \
//  std::cout << c << std::endl;              \
//  for (int j = 0; j < d_aij[i].size(); ++j) \
//    std::cout << d_aij[i][j] << std::endl;
#define db_print_row(i, c) void(0);

//----------------------------------------------------------------------------//
void Matrix::assemble()
{

  /*
   *  We construct in COO format, ie with (i,j,v) triplets.  This makes
   *  adding entries much easier.  Now, we need to order everything
   *  so that the resulting CSR storage has for all i-->j, with j in
   *  order.  The diagonal pointer is also set.  If (i,i,v) doesn't exist,
   *  we insert it, since that's needed for Gauss-Seidel, preconditioning,
   *  and other things.
   */

  Insist(!d_is_ready, "This matrix must not have been assembled already.");
  Insist(d_allocated, "This matrix must be allocated before assembling.");

  typedef std::vector<triplet_T>::iterator it_T;

  // sort the coo structure
  d_nnz = 0;
  for (int i = 0; i < d_m; i++)
  {
    db_print_row(i, "before sort")

    // do the sort
    std::sort(d_aij[i].begin(), d_aij[i].end(), compare_triplet);

    db_print_row(i, "after sort")

    // remove empty entries (if not entered, i,j==-1)
    while ((d_aij[i].end()-1)->j == -1)
      d_aij[i].pop_back();

    db_print_row(i, "after pop")

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
      {
        d_aij[i].insert(d_aij[i].begin()+d, triplet_T(i, i, 0.0));
        std::sort(d_aij[i].begin(), d_aij[i].end(), compare_triplet);

      }
    }

    db_print_row(i, "after diag")

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

//----------------------------------------------------------------------------//
void Matrix::clear()
{
  // save the memory used and reallocate

//  if (!d_allocated) return;
//  d_allocated = false;

  if (!d_is_ready) return;

  d_allocated = false;

  std::vector<int> nnz(d_m, 0);
  for (int i = 0; i < d_m; ++i)
    nnz[i] =  d_rows[i + 1] - d_rows[i];
  d_counter.clear();
  preallocate(&nnz[0]);

  if (!d_is_ready) return;

#ifdef CALLOW_ENABLE_PETSC
  int ierr = MatDestroy(&d_petsc_matrix);
  Assert(ierr == 0);
#endif

  delete [] d_values;
  delete [] d_columns;
  delete [] d_rows;
  delete [] d_diagonals;
  d_is_ready = false;
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of Matrix.cc
//----------------------------------------------------------------------------//
