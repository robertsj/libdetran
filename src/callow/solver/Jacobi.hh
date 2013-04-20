//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Jacobi.hh
 *  @brief  Jacobi class definition
 *  @author Jeremy Roberts
 *  date   Sep 13, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_JACOBI_HH_
#define callow_JACOBI_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  @class Jacobi
 *  @brief Uses Jacobi iteration to solve a system
 *
 *  Jacobi iteration in matrix form uses a splitting of
 *  the form
 *  @f[
 *      \mathbf{A} = \mathbf{L} + \mathbf{U} +  \mathbf{D} \, ,
 *  @f]
 *  which are strictly lower and upper triangle and diagonal,
 *  respectively.  The Jacobi iteration is then
 *  @f[
 *      \mathbf{D} x^{n+1} = -(\mathbf{L}+\mathbf{U})x^{n} + b
 *  @f]
 *  or
 *  @f[
 *      x^{n+1} = \overbrace{-\mathbf{D}^{-1}(\mathbf{L}+\mathbf{U})}^
 *                          {\mathbf{M}}x^{n} + \mathbf{D}^{-1}b \, .
 *  @f]
 *  The procedure converges if the iteration matrix is bounded, i.e.
 *  @f[
 *     \rho(\mathbf{M}) < 1 \, ,
 *  @f]
 *  which is guaranteed if \f$ \mathbf{A} \f$ is strictly diagonally
 *  dominant, meaning that
 *  @f[
 *      |a_{ii}| > \sum_{j \neq i} |a_{ij}| \, .
 *  @f]
 *
 *  Because we use a sparse matrix, we actually access the elements
 *  directly rather than via indexing.
 *
 */
class Jacobi: public LinearSolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Jacobi(const double atol,
         const double rtol,
         const int maxit,
         bool successive_norm = false);

  virtual ~Jacobi(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver::d_absolute_tolerance;
  using LinearSolver::d_relative_tolerance;
  using LinearSolver::d_maximum_iterations;
  using LinearSolver::d_residual;
  using LinearSolver::d_number_iterations;
  using LinearSolver::d_A;
  using LinearSolver::d_P;
  using LinearSolver::d_norm_type;

  // Use norm of residual (false=default) or successive iterates
  bool d_successive_norm;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  \param b  right hand side
   *  \param x  unknown vector
   */
  void solve_impl(const Vector &b, Vector &x);

};

} // end namespace callow

// Inline member definitions
#include "Jacobi.i.hh"

#endif // callow_JACOBI_HH_

//---------------------------------------------------------------------------//
//              end of file Jacobi.hh
//---------------------------------------------------------------------------//
