//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCILUT.hh
 *  @brief PCILUT class definition
 */
//----------------------------------------------------------------------------//

#ifndef callow_PCILUT_HH_
#define callow_PCILUT_HH_

#include "Preconditioner.hh"
#include "callow/matrix/Matrix.hh"

namespace callow
{

/**
 *  @class PCILUT
 *  @brief Implements the ILUT(p,t) preconditioner.
 *
 *  Following Saad, Algorithm 3, the ILUT(p, t) factorization is
 *
 *  @code
 *    given A, p, and t
 *    initialize L = eye(size(A)) and U = A
 *    row[0:n] = 0
 *    for i = 1:n
 *      row[0:n] = sparsecopy(a[i, :])
 *      for k = 0:i-1
 *        continue if row[k] == 0
 *        row[k] /= U[k, k]
 *        row[k] = drop_t(row[k], t)
 *        continue if row[k] == 0
 *        for j = k+1:n
 *           row[j] = sparseadd(row[j], -row[k]*U[k, j])
 *        end j
 *      end k
 *      row[0:n] = drop_p(row[0:n], p, t)
 *      L[i, 0:i-1] = sparsecopy(row[0:i-1])
 *      U[i, i:n] = sparsecopy(row[i:n])
 *      row[0:n] = 0
 *    end i
 *  @endcode
 *
 * Here, the `drop_t(v, t)` function returns 0 for `abs(v) < t` while
 * `drop_p(v[:], p, t)` limits the number of nonzero values in the row.
 * Various ways to implement the latter function are possible.  Saad
 * suggests first eliminating for the t threshold and then keeping the
 * largest p values in the separate L and U parts of the row.
 *
 * Reference:
 *   Saad, Y. "ILUT: A dual threshold incomplete LU factorization."
 *       Numerical Linear Algebra with Applications.  1.4 387-402 (1994)
 * (though this preprint was used:
 *  https://www-users.cs.umn.edu/~saad/PDF/umsi-92-38.pdf)
 */

class CALLOW_EXPORT PCILUT: public Preconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Preconditioner                    Base;
  typedef Base::SP_preconditioner           SP_preconditioner;
  typedef MatrixBase::SP_matrix             SP_matrix;
  typedef Matrix::SP_matrix                 SP_matrixfull;
  typedef Vector::SP_vector                 SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Construct an ILUT preconditioner for the explicit matrix A
  PCILUT(SP_matrix A, const size_t p = 1000, const double t = 0);

  /// SP constructor
  static SP_preconditioner
  Create(SP_matrix A, const size_t p = 1000, const double t = 0);

  /// Virtual destructor
  virtual ~PCILUT(){};

  /// Print my matrix.
  virtual void display(const std::string &name);

  /// Get my matrix.
  SP_matrixfull matrix() {return d_P;}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// Solve Px = b
  void apply(Vector &b, Vector &x);

private:

  /// Number of levels on either side of diagonal
  size_t d_p;
  /// Threshold for keeping an element of ILU
  double d_t;
  /// ILU decomposition of A
  SP_matrixfull d_P;
  /// Working vector
  Vector d_y;

};

} // end namespace callow

#endif // callow_PCILUT_HH_

//----------------------------------------------------------------------------//
//              end of file PCILUT.hh
//----------------------------------------------------------------------------//
