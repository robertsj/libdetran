//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OperatorShell.hh
 * \brief  OperatorShell class definition
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef OPERATORSHELL_HH_
#define OPERATORSHELL_HH_

#include "Operator.hh"
#include <iostream>

namespace detran
{

// Forward declare the wrapper
PetscErrorCode shell_multiply_wrapper(Mat A, Vec x, Vec y);

/*!
 *  \class OperatorMatrix
 *  \brief An operator with a matrix-free representation
 *
 *  All clients must implement the shell operators.
 */
class OperatorShell: public Operator
{

public:

  OperatorShell(const size_t m,
                const size_t n,
                void* context);


  /// Assemble the matrix.
  void assemble()
  {
    // Shells have nothing to assemble by default
  }

  /*!
   *  \brief Display the matrix to screen (or to output)
   *  \param  output    Flag indicating (stdout=0, ascii=1, binary=2)
   *  \param  name      File name for ascii or binary file
   */
  virtual void display(const int output = 0,
                       const std::string name) const
  {
    std::cout << "OperatorShell:" << std::endl
              << "  # rows = " << d_number_rows << std::endl
              << "  # cols = " << d_number_columns << std::endl
              << std::endl;
  }

protected:

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  /*!
   *  PETSc calls this function when applying the operator.  The function
   *  then calls the shell multiply function, when must be implemented
   *  by the derived class.
   */
  friend PetscErrorCode shell_multiply_wrapper(Mat A, Vec x, Vec y);

  /*!
   *  \brief Matrix-vector multiplication
   *  \param x  Input vector
   *  \param y  Output vector
   */
  virtual PetscErrorCode shell_multiply(Vec x, Vec y) = 0;

  /*!
   *  \brief Matrix-vector multiplication using matrix transpose.
   *  \param x  Input vector
   *  \param y  Output vector
   */
  virtual PetscErrorCode shell_multiply_transpose(Vec x, Vec y) = 0;


  /*!
   *  This must be called during construction of the base class.  It
   *  gives PETSc the wrapper functions.
   */
  PetscErrorCode set_operation();

};

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTION AND SETTER
//---------------------------------------------------------------------------//

PetscErrorCode shell_multiply_wrapper(Mat A, Vec x, Vec y)
{

  // Get the context and cast as MatrixShell
  PetscErrorCode ierr;
  void *context;
  ierr = MatShellGetContext(A, &context);
  Assert(!ierr);

  OperatorShell *foo = (OperatorShell*) context;
  // Call the actual apply operator.
  return foo->shell_multiply(x, y);
}

PetscErrorCode OperatorShell::set_operation()
{
  return MatShellSetOperation(d_A, MATOP_MULT,
                             (void(*)(void)) shell_multiply_wrapper);
}

} // end namespace detran

#endif // OPERATORSHELL_HH_ 

//---------------------------------------------------------------------------//
//              end of file OperatorShell.hh
//---------------------------------------------------------------------------//
