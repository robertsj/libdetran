//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussSeidel.hh
 *  @author robertsj
 *  @date   Sep 14, 2012
 *  @brief  GaussSeidel class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_GAUSSSEIDEL_HH_
#define callow_GAUSSSEIDEL_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  @class GaussSeidel
 *  @brief Uses Gauss-Seidel iteration to solve a system
 *
 *  Gauss-Seidel iteration in matrix form uses a splitting of
 *  the form
 *  @f[
 *      \mathbf{A} = \mathbf{L} + \mathbf{U} +  \mathbf{D} \, ,
 *  @f]
 *  which are strictly lower and upper triangle and diagonal,
 *  respectively.  The Gauss-Seidel iteration is then
 *  @f[
 *      \mathbf{D + L} x^{n+1} = -\mathbf{U}x^{n} + b
 *  @f]
 *  or
 *  @f[
 *      x^{n+1} = \overbrace{-\mathbf{D+L}^{-1}(\mathbf{U})}^
 *                          {\mathbf{M}}x^{n} + \mathbf{D+L}^{-1}b \, .
 *  @f]
 *  Alternatively, one can swap @f$ U @f$ and @f$ L @f$ to produce
 *  the backward Gauss-Seidel iteration.  If used together, one
 *  has symmetric Gauss-Seidel iteration.
 *
 *  It can be shown that Gauss-Seidel converges if \ref Jacobi
 *  converges.
 *
 *  Because we use a sparse matrix, we actually access the elements
 *  directly rather than via indexing.
 *
 */
template<class T>
class GaussSeidel: public LinearSolver<T>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver<T> Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  GaussSeidel(const double atol, const double rtol, const int maxit);

  virtual ~GaussSeidel(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver<T>::d_absolute_tolerance;
  using LinearSolver<T>::d_relative_tolerance;
  using LinearSolver<T>::d_maximum_iterations;
  using LinearSolver<T>::d_L1_residual;
  using LinearSolver<T>::d_L2_residual;
  using LinearSolver<T>::d_LI_residual;
  using LinearSolver<T>::d_number_iterations;
  using LinearSolver<T>::d_A;
  using LinearSolver<T>::d_P;
  using LinearSolver<T>::d_norm_type;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  @param b  right hand side
   *  @param x  unknown vector
   */
  void solve_impl(const Vector<T> &b, Vector<T> &x);

};

} // end namespace callow

// Inline member definitions
#include "GaussSeidel.i.hh"


#endif /* callow_GAUSSSEIDEL_HH_ */
