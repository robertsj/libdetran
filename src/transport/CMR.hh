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
 *  Consider the one group problem
 *  \f[
 *      \Big ( \mu  \frac{\partial}{\partial x}
 *           + \eta \frac{\partial}{\partial y}
 *           + \xi  \frac{\partial}{\partial z}
 *      \Big ) \psi + \Sigma_t(\vec{r}) \psi(\vec{r}, \hat{\Omega}) =
 *      \Sigma_s(\vec{r}) \phi(\vec{r}) + q(\vec{r}, \hat{\Omega}) \, .
 *  \f]
 *  where we have explicitly separated the within-group scatter
 *  from the external source.
 *
 *  Suppose we integrate this over a coarse cell of volume
 *  \f$ V = \Delta_x \Delta_y \Delta_z \f$ and over the
 *  angular space, \f$ 4\pi \f$. The
 *  streaming terms give rise to terms like
 *
 *  \f[
 *  \begin{split}
 *      & \int_{4\pi} \mu \int^{\Delta_y}_{0} \int^{\Delta_z}_{0} dy dz
 *         \Bigg ( \frac{\partial \psi}{\partial x} \Big |_{x=\Delta_x}
 *                -\frac{\partial \psi}{\partial x} \Big |_{x=0} \Bigg ) =
 *      \int^{\Delta_y}_{0} \int^{\Delta_z}_{0} dy dz
 *        \Big ( J_x(\Delta_x, y, z) - J_x(0, y, z) \Big ) \, ,
 *  \end{split}
 *  \f]
 *
 *  where \f$ J_x \f$ is the \f$x\f$-directed partial current.  Similar
 *  terms can be found for the other directions.
 *
 *  We can write the integrated transport equation as
 *
 *  \f[
 *     \Delta_y \Delta_z \Big ( \bar{J}_x(\Delta_x) - \bar{J}_x(0) \Big ) +
 *     \Delta_x \Delta_z \Big ( \bar{J}_y(\Delta_y) - \bar{J}_y(0) \Big ) +
 *     \Delta_x \Delta_y \Big ( \bar{J}_z(\Delta_z) - \bar{J}_z(0) \Big ) +
 *     \Delta_x \Delta_y \Delta_z \bar{\Sigma_r} \bar{\phi}(\vec{r}) =
 *     \Delta_x \Delta_y \Delta_z \bar{q}(\vec{r})
 *  \f]
 *
 *  where we have  defined the average current across a face of the box as
 *
 *  \f[
 *     \bar{J}_x(x) =
 *       \frac{\int^{\Delta_y}_{0} \int^{\Delta_z}_{0} dy dz J_x(x, y, z)}
 *            {\int^{\Delta_y}_{0} \int^{\Delta_z}_{0} dy dz} \, .
 *  \f]
 *
 *  and similar averages for the removal rate and source.  Note, it is
 *  the \e removal rate of interest on the left.  This requires a coarse
 *  mesh averaged removal cross section and scalar flux.  The removal
 *  cross section is simply the total cross section minus the
 *  within-group scattering cross section.  As a consequence, the
 *  source \f$ q \f$ does \e not include the within-group scatter source.
 *
 *  We simplify notation by writing
 *
 *  \f[
 *     \sum_{p \in (x,y,z)} \frac{1}{\Delta_{p}}
 *       \Big ( \bar{J}_{p}(\Delta_{p}) - \bar{J}_{p}(0) \Big ) +
 *     \bar{\Sigma_r} \bar{\phi}(\vec{r}) =
 *     \bar{q}(\vec{r})
 *  \f]
 *
 *  This neutron balance equation is not satisfied when the flux
 *  in not converged.  However, if we can somehow force the flux
 *  to satisfy balance at the coarse level, then we can apply
 *  any correction to the fine level.  Let us define rebalance
 *  factors \f$ f_i \f$ such that
 *
 *  \f[
 *     \sum_{p \in (x,y,z)} \frac{1}{\Delta_{p}}
 *       \Big ( \bar{J}_{p}(\Delta_{p}) - \bar{J}_{p}(0) \Big ) +
 *     f_i \bar{\Sigma_r} \bar{\phi}_i =
 *     \bar{q}(\vec{r})
 *  \f]
 *
 */
template <class D>
class CMR : public WithinGroupAcceleration<D>
{

public:

  /// \name Useful Typedefs
  // \{

  typedef SP<CMR>                         SP_acceleration;
  typedef WithingroupAcceleration<D>      Base;
  typedef typename Base::SP_acceleration  SP_base;
  typedef Base::SP_input                  SP_input;
  typedef Base::SP_material               SP_material;
  typedef Base::SP_coarsemesh             SP_coarsemesh;
  typedef Base::SP_currenttally           SP_currenttally;
  typedef typename Base::SP_sweepsource   SP_sweepsource;

  // \}

  /*!
   *  \brief Constructor
   *
   *  \param input          Input database
   *  \param material       Material database
   *  \param coarsemesh     Coarse mesh
   *  \param currenttally   Current tally
   */
  CMR(SP_input input,
      SP_material material,
      SP_coarsemesh coarsemesh,
      SP_currenttally currenttally);

  /// \name Public Interface
  /// \{


  /*!
   *  \brief Solve the low order equation to update the scalar flux moments
   *
   *  \param  group   Energy group for this solve
   *  \param  phi     Reference to the group flux to be updated.
   *  \param  source  Smart pointer to up-to-date sweep source.
   *
   */
  void accelerate(u_int group,
                  State::moments_type &phi,
                  SP_sweepsource source);

  /// \}

private:

  /// \name Private Data
  /// \{

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
  void integrate(u_int group,
                 State::moments_type &phi,
                 SP_sweepsource source);

  /// \}

};


} // end namespace detran

#include "CMR.i.hh"

#endif /* CMR_HH_ */
