//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PowerIteration.hh
 *  @brief  PowerIteration
 *  @author Jeremy Roberts
 *  @date   Sep 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef POWERITERATION_HH_
#define POWERITERATION_HH_

#include "EigenSolver.hh"

namespace callow
{

/**
 *  @class PowerIteration
 *  @brief Solve the eigenvalue problem with the power method
 */
template <class T>
class PowerIteration: public EigenSolver<T>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef EigenSolver<T>                    Base;
  typedef typename Base::SP_matrix          SP_matrix;
  typedef typename Base::SP_solver          SP_solver;
  typedef typename Base::SP_vector          SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  PowerIteration(const double    tol = 1e-6,
                 const int       maxit = 100);

  virtual ~PowerIteration(){}

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base members
  using Base::d_tolerance;
  using Base::d_maximum_iterations;
  using Base::d_name;
  using Base::d_residual_norm;
  using Base::d_number_iterations;
  using Base::d_A;
  using Base::d_B;
  using Base::d_solver;
  using Base::d_monitor_level;
  using Base::d_lambda;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  virtual void solve_impl(Vector<T> &x, Vector<T> &x0);

};

} // end namespace callow

#include "PowerIteration.i.hh"

#endif // POWERITERATION_HH_ 

//---------------------------------------------------------------------------//
//              end of file PowerIteration.hh
//---------------------------------------------------------------------------//
