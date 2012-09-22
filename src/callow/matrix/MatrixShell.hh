//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixShell.hh
 * \author robertsj
 * \date   Sep 18, 2012
 * \brief  MatrixShell class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXSHELL_HH_
#define callow_MATRIXSHELL_HH_

#include "MatrixBase.hh"

namespace callow
{

/*!
 *  \class MatrixShell
 *  \brief Defines a matrix free operator
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
template <class T>
class MatrixShell: public MatrixBase<T>
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<MatrixShell<T> >    SP_matrix;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  /// the context is the "this" of the caller
  MatrixShell(void* context);
  /// construct with known sizes
  MatrixShell(void* context, const int m, const int n);
  virtual ~MatrixShell();

  void set_size(const int m, const int n)
  {
    Require(m > 0 and n > 0);
    d_m = m;
    d_n = n;
    d_sizes_set = true;
#ifdef DETRAN_ENABLE_PETSC
    PetscErrorCode ierr;
    ierr = MatCreate(PETSC_COMM_SELF, &d_petsc_matrix);
    ierr = MatSetSizes(d_petsc_matrix, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
    ierr = MatSetType(d_petsc_matrix, MATSHELL);
    ierr = MatShellSetContext(d_petsc_matrix, d_context);
    Ensure(!ierr);
#endif
    set_operation();
    d_is_ready = true;
  }

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT
  //---------------------------------------------------------------------------//

  // default shell assemble does nothing
  void assemble()
  {
    /* ... */
  }
  // default shell display gives just the sizes
  void display() const
  {
    std::cout << "MatrixShell:" << std::endl
              << "  # rows = " << d_m << std::endl
              << "  # cols = " << d_n << std::endl
              << std::endl;
  }
  // the client must implement the action y <-- A * x
  virtual void multiply(const Vector<T> &x,  Vector<T> &y) = 0;
  // the client must implement the action y <-- A' * x
  virtual void multiply_transpose(const Vector<T> &x, Vector<T> &y) = 0;

protected:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  /// expose base members
  using MatrixBase<T>::d_m;
  using MatrixBase<T>::d_n;
  using MatrixBase<T>::d_sizes_set;
  using MatrixBase<T>::d_is_ready;

#ifdef DETRAN_ENABLE_PETSC
  using MatrixBase<T>::d_petsc_matrix;
#endif

private:

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

} // end namespace callow

#include "MatrixShell.i.hh"

#endif /* callow_MATRIXSHELL_HH_ */
