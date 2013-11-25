//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixDense.hh
 *  @brief MatrixDense class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_MATRIXDENSE_HH_
#define callow_MATRIXDENSE_HH_

#include "MatrixBase.hh"
#include <ostream>
#include <vector>
#include <string>

namespace callow
{

#define MATRIXDENSE_COLMAJ

/**
 *  @class MatrixDense
 *  @brief Dense matrix
 *
 *  This class uses a column-major storage format lie Fortran.  This makes
 *  it easier to couple with PETSc, BLAS, etc.
 */

class CALLOW_EXPORT MatrixDense: public MatrixBase
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum insert_type
  {
    INSERT, ADD, END_INSERT_TYPE
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<MatrixDense>  SP_matrix;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  // construction with sizing but deferred allocation
  MatrixDense(const int m, const int n, const double v = 0.0);
  // copy constructor
  MatrixDense(const MatrixDense &A);
  // destructor
  virtual ~MatrixDense();
  // sp constructor
  static SP_matrix Create(const int m, const int n, const double v = 0.0);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// add one value (return false if can't add)
  bool insert(int  i, int  j, double  v, const int type = INSERT);
  /// add a row  (return false if can't add)
  bool insert_row(int  i, double *v, const int type = INSERT);
  /// add a column  (return false if can't add)
  bool insert_col(int  j, double *v, const int type = INSERT);

  /// value at a cardinal index
  const double& operator[](const int p) const;
  double& operator[](const int p);

  /// value at ij and returns 0 if not present
  const double& operator()(const int i, const int j) const;
  double& operator()(const int i, const int j);

  // get underlying storage and indexing. careful!
  double* values() {return d_values;}

  /// print (i, j, v) to ascii file with 1-based indexing for matlab
  void print_matlab(std::string filename = "matrix.out") const;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  // doesn't do anything for dense matrix
  void assemble(){ /* ... */ }
  // action y <-- A * x
  void multiply(const Vector &x,  Vector &y);
  // action y <-- A' * x
  void multiply_transpose(const Vector &x, Vector &y);
  // pretty print to screen
  void display(bool forceprint = false) const;
  // clear the matrix contents
  void clear();

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// expose base members
  using MatrixBase::d_m;
  using MatrixBase::d_n;
  using MatrixBase::d_is_ready;

#ifdef CALLOW_ENABLE_PETSC
  using MatrixBase::d_petsc_matrix;
#endif

  /// matrix elements
  double* d_values;

private:

};

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<MatrixDense>)

} // end namespace callow

// Inline members
#include "MatrixDense.i.hh"

#endif // callow_MATRIXDENSE_HH_

//----------------------------------------------------------------------------//
//              end of file MatrixDense.hh
//----------------------------------------------------------------------------//
