//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIteration.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  PowerIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef POWERITERATION_HH_
#define POWERITERATION_HH_

// Detran
#include "Eigensolver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class PowerIteration
 * \brief Solves the eigenvalue problem via the power method.
 *
 * The eigenvalue problem can be cast in the form
 * \f[
 *     \mathbf{A}d = kd \,
 * \f]
 * where \f$ d \f$ is the fission density and \f$ k \f$ is the
 * eigenvalue.  See Eigensolver for more details on this formulation.
 *
 * The power method solves the eigenproblem using the iteration
 * \f[
 *     d^{l+1} \leftarrow \mathbf{A} d^{l} / k^{l}
 * \f]
 * and
 * \f[
 *     k^{l+1} \leftarrow || d^{l+1} || \, .
 * \f]
 * Traditionally, the norm used to define the updated eigenvalue
 * is a fission-weighted sum of the group fluxes, which is
 * equivalent to an L1 norm of the density if all the fluxes
 * and fission cross sections are positive, true for all
 * physical problems (and barring numerical issues due to
 * discretization).  Here, we use the L1 norm, and begin
 * with \f$ || d^{0} || = 1 \f$.
 *
 */
//---------------------------------------------------------------------------//

template <class D>
class PowerIteration: public Eigensolver<D>
{

public:

  typedef SP<PowerIteration<D> >                SP_solver;
  typedef Eigensolver<D>                        Base;
  typedef typename Base::SP_solver              SP_base;
  typedef typename
      MultigroupSolver<D>::SP_solver            SP_mg_solver;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  typedef FissionSource::SP_source              SP_fissionsource;

  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param fission_source    Fission source.
   */
  PowerIteration(SP_input           input,
                 SP_state           state,
                 SP_mesh            mesh,
                 SP_material        material,
                 SP_quadrature      quadrature,
                 SP_boundary        boundary,
                 SP_fissionsource   q_f);

  /// SP Constructor
  static SP<PowerIteration<D> >
  Create(SP<detran::InputDB>          input,
         SP<detran::State>            state,
         SP<detran::Mesh>             mesh,
         SP<detran::Material>         material,
         SP<detran::Quadrature>       quadrature,
         SP<detran::BoundaryBase<D> > boundary,
         SP<detran::FissionSource>    q_f)
  {
    SP_solver p;
    p = new PowerIteration(input, state, mesh, material,
                           quadrature, boundary, q_f);
    return p;
  }

  /// Solve the eigenvalue problem.
  void solve();

  /// Unimplemented DBC function.
  bool is_valid() const
  {
    return true;
  }


protected:

  using Base::b_input;
  using Base::b_state;
  using Base::b_mesh;
  using Base::b_material;
  using Base::b_quadrature;
  using Base::b_boundary;
  using Base::b_fissionsource;
  using Base::b_mg_solver;
  using Base::b_max_iters;
  using Base::b_tolerance;
  using Base::b_print_out;
  using Base::b_print_interval;

  /// \name Protected Data
  /// \{

  /// Display Aitken extrapolation
  bool d_aitken;

  /// \}

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "PowerIteration.i.hh"

#endif /* POWERITERATION_HH_ */
