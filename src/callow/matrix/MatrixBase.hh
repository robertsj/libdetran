//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MatrixBase.hh
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  MatrixBase class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXBASE_HH_
#define callow_MATRIXBASE_HH_

#include "utilities/Definitions.hh"
#include "callow/callow_config.hh"
#include "callow/vector/Vector.hh"
#include "utilities/SP.hh"

namespace callow
{

class CALLOW_EXPORT MatrixBase
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<MatrixBase >     SP_matrix;
  typedef Vector::SP_vector                     SP_vector;
  typedef detran_utilities::size_t              size_t;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  MatrixBase()
    : d_m(0)
    , d_n(0)
    , d_sizes_set(false)
    , d_is_ready(false)
  {
    /* ... */
  }

  MatrixBase(const int m, const int n)
    : d_m(m)
    , d_n(n)
    , d_is_ready(false)
  {
    Require(m > 0 && n > 0);
    d_sizes_set = true;
  }

  virtual ~MatrixBase()
  {
#ifdef CALLOW_ENABLE_PETSC
    // destroy the petsc matrix.  note, since we constructed it
    // using our own arrays, those still need to be deleted.
    MatDestroy(&d_petsc_matrix);
#endif
  }

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  virtual void set_size(const int m, const int n)
  {
    Require(m > 0 && n > 0);
    d_m = m;
    d_n = n;
    d_sizes_set = true;
  }

  int number_rows() const { return d_m; }
  int number_columns() const { return d_n; }
  bool is_ready() const { return d_is_ready; }

#ifdef CALLOW_ENABLE_PETSC
  Mat petsc_matrix() {return d_petsc_matrix;}
#endif

  // multiply with SP vectors
  void multiply(SP_vector x,  SP_vector y)
  {
    multiply(*x, *y);
  }
  // multiply transpose with SP vectors
  void multiply_transpose(SP_vector x, SP_vector y)
  {
    multiply_transpose(*x, *y);
  }

  // compute and print the explicit operator (even if shell)
  virtual void compute_explicit(std::string filename = "matrix.out")
  {
#ifdef CALLOW_ENABLE_PETSC
  Mat A;
  MatComputeExplicitOperator(d_petsc_matrix, &A);
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_WRITE, &viewer);
  MatView(A, viewer);
  PetscViewerDestroy(&viewer);
#endif
  }

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT
  //---------------------------------------------------------------------------//

  // postprocess storage
  virtual void assemble() = 0;

  // action y <-- A * x
  virtual void multiply(const Vector &x,  Vector &y) = 0;
  // action y <-- A' * x
  virtual void multiply_transpose(const Vector &x, Vector &y) = 0;

  // pretty print to screen
  virtual void display(bool forceprint = false) const = 0;
  // print output for reading into matlab
  virtual void print_matlab(std::string filename = "matrix.out") const
  {
    /* ... */
  }

protected:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  /// number of rows
  int d_m;
  /// number of columns
  int d_n;
  /// are m and n set?
  bool d_sizes_set;
  /// am i good to go?
  bool d_is_ready;

#ifdef CALLOW_ENABLE_PETSC
  Mat d_petsc_matrix;
#endif

};

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<MatrixBase>)

} // end namespace callow

#endif /* callow_MATRIXBASE_HH_ */
