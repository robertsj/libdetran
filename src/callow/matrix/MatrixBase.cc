//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixBase.cc
 *  @brief MatrixBase member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MatrixBase.hh"
#include "MatrixDense.hh"
#include "Vector.hh"
#include "utils/Typedefs.hh"
#include <iostream>

namespace callow
{

//----------------------------------------------------------------------------//
MatrixBase::MatrixBase()
  : d_m(0)
  , d_n(0)
  , d_sizes_set(false)
  , d_is_ready(false)
  , d_petsc_matrix(NULL)
{
  /* ... */
}

//----------------------------------------------------------------------------//
MatrixBase::MatrixBase(const int m, const int n)
  : d_m(m)
  , d_n(n)
  , d_is_ready(false)
  , d_petsc_matrix(NULL)
{
  Require(m > 0 && n > 0);
  d_sizes_set = true;
}

//----------------------------------------------------------------------------//
MatrixBase::~MatrixBase()
{
#ifdef CALLOW_ENABLE_PETSC
  // destroy the petsc matrix.  note, since we constructed it
  // using our own arrays, those still need to be deleted.
  MatDestroy(&d_petsc_matrix);
#endif
}

//----------------------------------------------------------------------------//
void MatrixBase::set_size(const int m, const int n)
{
  Require(m > 0 && n > 0);
  d_m = m;
  d_n = n;
  d_sizes_set = true;
}

//----------------------------------------------------------------------------//
void MatrixBase::compute_explicit(std::string filename)
{
#ifdef CALLOW_ENABLE_PETSC
Mat A;
MatComputeExplicitOperator(d_petsc_matrix, &A);
PetscViewer viewer;
PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_WRITE,
                      &viewer);
MatView(A, viewer);
PetscViewerDestroy(&viewer);
#endif
  if (true)
  {
	// Allocate dense matrix for explicit construction
	MatrixDense A(d_m, d_n, 0.0);
    for (int j = 0; j < d_n; ++j)
    {
    	// Define e_i w
    	Vector e(d_n, 0.0);
    	e[j] = 1.0;

    	// Define output
    	Vector a_j(d_m, 0.0);

    	// a = A*e
        this->multiply(e, a_j);

        // Set the column
        for (int i = 0; i < d_m; ++i)
        {
          A(i, j) = a_j[i];
        }
    }
    // Output the dense matrix
    A.print_matlab(filename);
  }
}

//----------------------------------------------------------------------------//
void MatrixBase::print_matlab(std::string filename) const
{
  /* ... */
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of MatrixBase.cc
//----------------------------------------------------------------------------//
