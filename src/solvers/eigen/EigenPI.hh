//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenPI.hh
 *  @brief EigenPI class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENPI_HH_
#define detran_EIGENPI_HH_

#include "Eigensolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class EigenPI
 *  @brief Solves the eigenvalue problem via the power method.
 *
 *  The eigenvalue problem can be cast in the form
 *  @f[
 *      \mathbf{A}d = kd \,
 *  @f]
 *  where \f$ d \f$ is the fission density and \f$ k \f$ is the
 *  eigenvalue.  See \ref Eigensolver for more details on this formulation.
 *
 *  The power method solves the eigenproblem using the iteration
 *  @f[
 *      d^{l+1} \leftarrow \mathbf{A} d^{l} / k^{l}
 *  @f]
 *  and
 *  @f[
 *      k^{l+1} \leftarrow || d^{l+1} || \, .
 *  @f]
 *  Traditionally, the norm used to define the updated eigenvalue
 *  is a fission-weighted sum of the group fluxes, which is
 *  equivalent to an L1 norm of the density if all the fluxes
 *  and fission cross sections are positive, true for all
 *  physical problems (and barring numerical issues due to
 *  discretization).  Here, we use the L1 norm, and begin
 *  with \f$ || d^{0} || = 1 \f$.
 *
 *  Note, this is a hand-coded power iteration implementation that
 *  can be used with nonlinear acceleration.
 *
 */
//----------------------------------------------------------------------------//

template <class D>
class EigenPI: public Eigensolver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Eigensolver<D>                            Base;
  typedef typename Base::SP_solver                  SP_solver;
  typedef typename Base::SP_mg_solver               SP_mg_solver;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_fissionsource           SP_fissionsource;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver         Multigroup solver
   */
  EigenPI(SP_mg_solver mg_solver);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the eigenvalue problem.
  void solve();

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Expose base members
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_fissionsource;
  using Base::d_number_groups;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;
  using Base::d_mg_solver;

  /// Display Aitken extrapolation
  bool d_aitken;

  /// Over-relaxation parameter
  double d_omega;

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "EigenPI.i.hh"

#endif /* detran_EIGENPI_HH_ */

//----------------------------------------------------------------------------//
//              end of file EigenPI.hh
//----------------------------------------------------------------------------//
