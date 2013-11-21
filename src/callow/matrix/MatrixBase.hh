//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixBase.hh
 *  @brief MatrixBase class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

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

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<MatrixBase >     SP_matrix;
  typedef Vector::SP_vector                     SP_vector;
  typedef detran_utilities::size_t              size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// default constructor
  MatrixBase();
  /// constructor with explicit sizing
  MatrixBase(const int m, const int n);
  /// virtual destructor
  virtual ~MatrixBase();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// set the matrix sizes
  virtual void set_size(const int m, const int n);
  /// get matrix sizes
  //@{
  int number_rows() const { return d_m; }
  int number_columns() const { return d_n; }
  //@}
  /// is the matrix ready to use?
  bool is_ready() const { return d_is_ready; }
  /// return a PETSc matrix
  Mat petsc_matrix() {return d_petsc_matrix;}

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
  virtual void compute_explicit(std::string filename = "matrix.out");
  // print output for reading into matlab
  virtual void print_matlab(std::string filename = "matrix.out") const;

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
  // clear the matrix contents
  virtual void clear() {}

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// number of rows
  int d_m;
  /// number of columns
  int d_n;
  /// are m and n set?
  bool d_sizes_set;
  /// am i good to go?
  bool d_is_ready;
  /// PETSc matrix
  Mat d_petsc_matrix;

};

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<MatrixBase>)

} // end namespace callow

#endif /* callow_MATRIXBASE_HH_ */

//----------------------------------------------------------------------------//
//              end of MatrixBase.hh
//----------------------------------------------------------------------------//
