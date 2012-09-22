//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Richardson.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Richardson class definition.
 */
//---------------------------------------------------------------------------//

#ifndef RICHARDSON_HH_
#define RICHARDSON_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  \class Richardson
 *  \brief Uses (modified) Richardson iteration to solve a system
 *
 *  Richardson iteration solves a linear system via the
 *  process
 *  \f[
 *     x^{(n+1)} = (\mathbf{I - \omega A})x^{(n)} + \omega b
 *  \f]
 *  where \f$ \omega \f$ is something like a relaxation factor
 *  that takes on values between (roughly) 0 and 2.  By default,
 *  \f$ \omega = 1 \f$.
 *
 */
template<class T>
class Richardson: public LinearSolver<T>
{

public:

  typedef LinearSolver<T> Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Richardson(const double atol, const double rtol, const int maxit,
             const double omega = 1.0);

  virtual ~Richardson(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  void set_omega(const double omega)
  {
    d_omega = omega;
  }

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

  /// relaxation factor
  double d_omega;

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
#include "Richardson.i.hh"

#endif /* RICHARDSON_HH_ */
