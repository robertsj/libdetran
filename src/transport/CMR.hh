/*
 * CMR.hh
 *
 *  Created on: May 17, 2012
 *      Author: robertsj
 */

#ifndef CMR_HH_
#define CMR_HH_

// Detran
#include "WithinGroupAcceleration.hh"

// System

namespace detran
{



/*!
 *  \class CMR
 *  \brief Coarse Mesh Rebalance
 *
 *  Describe me.
 *
 *  Note, this applies only to the within group problem.
 */
template <class D>
class CMR : public WithinGroupAcceleration<D>
{

public:

  /// \name Useful Typedefs
  // \{

  typedef SP<CMR>                                     SP_acceleration;
  typedef WithinGroupAcceleration<D>                  Base;
  typedef typename Base::SP_acceleration              SP_base;
  typedef typename Base::SP_mesh                      SP_mesh;
  typedef typename Base::SP_material                  SP_material;
  typedef typename Base::SP_quadrature                SP_quadrature;
  typedef typename Base::SP_state                     SP_state;
  typedef typename Base::face_flux_type               face_flux_type;
  typedef typename SweepSource<D>::SP_sweepsource     SP_sweepsource;

  // \}

  /*!
   *  \brief Constructor
   *
   *  \param mesh       Mesh smart pointer
   *  \param material   Material smart pointer
   *  \param quadrature Quadrature smart pointer
   */
  CMR(SP_mesh mesh, SP_material material, SP_quadrature quadrature);

  /// Virtual destructor.
  ~CMR(){}

  // PUBLIC INTERFACE

  /*!
   *  \brief Add contribution to an arbitrary function of the coarse
   *         mesh edge flux.
   *
   *  For CMR, we compute the coarse mesh edge partial currents.
   *
   *  \param  i   x mesh index
   *  \param  j   y mesh index
   *  \param  k   z mesh index
   *  \param  o   octant
   *  \param  a   angle within octant
   *  \param  psi edge angular flux
   */
  void tally(int i, int j, int k, int o, int a, face_flux_type psi);

  /*!
   *  \brief Solve the low order equation to update the scalar flux moments
   *
   *  The removal reaction rate and total source is integrated in each
   *  coarse mesh.  The partial currents at each coarse mesh surfaces
   *  are precomputed via the tally method during sweeps.  The system
   *  of equations for rebalance factors is produced and solved, leading
   *  to an updated fine mesh flux.
   *
   *  \param  phi     Reference to the group flux to be updated.
   *  \param  source  Smart pointer to up-to-date sweep source.
   *
   */
  void update(State::moments_type &phi, SP_sweepsource source);

  /*!
   *  \brief Create the coarse mesh and perform associated allocations.
   */
  void initialize(int level);

  void reset();


private:

  // Expose Base Class members
  using Base::b_mesh;
  using Base::b_coarse_mesh;

  /// \name Private Data
  /// \{

  /// Positive sense edge partial current.
  vec2_dbl d_J_pos;

  /// Negative sense edge partial current.
  vec2_dbl d_J_neg;

  /// Cell-integrated removal rate
  vec_dbl d_R;

  /// Cell-integrated source
  vec_dbl d_Q;

  /// Rebalance factors.
  vec2_int d_f;

  /// \}

  /// \name Implementation
  /// \{

  /*!
   *  \brief Compute the coarse mesh integrated removal rate and source
   *  \param  phi     Reference to fine mesh flux
   *  \param  source  Smart pointer to sweep source
   */
  void integrate(State::moments_type &phi, SP_sweepsource source);

  /// \}

};


} // end namespace detran

#include "CMR.i.hh"

#endif /* CMR_HH_ */
