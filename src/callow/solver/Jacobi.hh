//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Jacobi.hh
 * \brief  Jacobi class definition
 * \author Jeremy Roberts
 * \date   Sep 13, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_JACOBI_HH_
#define callow_JACOBI_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  \class Jacobi
 *  \brief Uses Jacobi iteration to solve a system
 *
 *  Jacobi iteration in matrix form uses a splitting of
 *  the form
 *  \f[
 *      \mathbf{A} = \mathbf{L} + \mathbf{U} +  \mathbf{D} \, ,
 *  \f]
 *  which are strictly lower and upper triangle and diagonal,
 *  respectively.  The Jacobi iteration is then
 *  \f[
 *      \mathbf{D} x^{n+1} = -(\mathbf{L}+\mathbf{U})x^{n} + b
 *  \f]
 *  or
 *  \f[
 *      x^{n+1} = \overbrace{-\mathbf{D}^{-1}(\mathbf{L}+\mathbf{U})}^
 *                          {\mathbf{M}}x^{n} + \mathbf{D}^{-1}b \, .
 *  \f]
 *  The procedure converges if the iteration matrix is bounded, i.e.
 *  \f[
 *     \rho(\mathbf{M}) < 1 \, ,
 *  \f]
 *  which is guaranteed if \f$ \mathbf{A} \f$ is strictly diagonally
 *  dominant, meaning that
 *  \f[
 *      |a_{ii}| > \sum_{j \neq i} |a_{ij}| \, .
 *  \f]
 *
 *  Because we use a sparse matrix, we actually access the elements
 *  directly rather than via indexing.
 *
 */
template<class T>
class Jacobi: public LinearSolver<T>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver<T> Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Jacobi(const double atol, const double rtol, const int maxit);

  virtual ~Jacobi(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver<T>::status;
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
   *  \param b  right hand side
   *  \param x  unknown vector
   */
  void solve_impl(const Vector<T> &b, Vector<T> &x);

};

} // end namespace callow

// Inline member definitions
#include "Jacobi.i.hh"

#endif // callow_JACOBI_HH_

//---------------------------------------------------------------------------//
//              end of file Jacobi.hh
//---------------------------------------------------------------------------//
