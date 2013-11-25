//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixShell.hh
 *  @brief MatrixShell class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_MATRIXSHELL_HH_
#define callow_MATRIXSHELL_HH_

#include "MatrixBase.hh"

namespace callow
{

/**
 *  @class MatrixShell
 *  @brief Defines a matrix free operator
 *
 *  For many iterative methods, only the action of the operator
 *  on a vector is required.  Frequently, constructing a matrix
 *  explicitly is too memory intensive, and so a matrix free
 *  operator that defines the action is desirable.
 *
 *  This class is abstract. The client must define the action of
 *  the matrix and its transpose (though for the latter, the
 *  client may simply throw an exception to forbid its use).
 *
 */
class CALLOW_EXPORT MatrixShell: public MatrixBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<MatrixShell>    SP_matrix;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// the context is the "this" of the caller
  MatrixShell(void* context);
  /// construct with known sizes
  MatrixShell(void* context, const int m, const int n);
  virtual ~MatrixShell();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// overloaded size setter so that the appropriate PETSc Mat is made
  void set_size(const int m, const int n = 0);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  // default shell assemble does nothing
  void assemble(){}
  // default shell display gives just the sizes
  void display(bool forceprint = false) const;
  // the client must implement the action y <-- A * x
  void multiply(const Vector &x,  Vector &y) = 0;
  // the client must implement the action y <-- A' * x
  void multiply_transpose(const Vector &x, Vector &y) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// context of my caller
  void* d_context;
  /// tell petsc about the operation, if applicable
  void set_operation();

};

#ifdef DETRAN_ENABLE_PETSC
// this is the function petsc actual calls; internally, it redirects
// to our own shell operation
inline PetscErrorCode shell_multiply_wrapper(Mat A, Vec x, Vec y);
#endif

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<MatrixShell>)

} // end namespace callow

#include "MatrixShell.i.hh"

#endif /* callow_MATRIXSHELL_HH_ */

//----------------------------------------------------------------------------//
//              end of MatrixShell.hh
//----------------------------------------------------------------------------//
