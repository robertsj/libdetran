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
  CMR(SP_input,
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
