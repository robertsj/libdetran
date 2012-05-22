/*
 * WithinGroupAcceleration.hh
 *
 *  Created on: May 18, 2012
 *      Author: robertsj
 */

#ifndef WITHINGROUPACCELERATION_HH_
#define WITHINGROUPACCELERATION_HH_

// Detran
#include "Acceleration.hh"
#include "SweepSource.hh"

namespace detran
{

/*!
 *  \class WithinGroupAcceleration
 *  \brief Base class for within-group coarse mesh acceleration schemes.
 *
 *
 */
template <class D>
class WithinGroupAcceleration : public Acceleration<D>
{

public:

  /// \name Useful Typedefs
  // \{

  typedef SP<WithinGroupAcceleration>                 SP_acceleration;
  typedef Acceleration<D>                             Base;
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
  WithinGroupAcceleration(SP_mesh mesh, SP_material material, SP_quadrature quadrature);

  /// Virtual destructor.
  virtual ~WithinGroupAcceleration(){};

  // \name Public Interface
  // \{

  /*!
   *  \brief Solve the low order equation to update the scalar flux moments
   *
   *  Each acceleration scheme has tallied the coarse mesh edge flux
   *  functions needed to perform the low order update.  The update
   *  function implements that update.  The fine mesh flux moments are
   *  passed, and the sweep source is passed to compute coarse mesh
   *  quantities.
   *
   *  \param  phi     Reference to the group flux to be updated.
   *  \param  source  Smart pointer to up-to-date sweep source.
   *
   */
  virtual void update(State::moments_type &phi, SP_sweepsource source) = 0;

};

} // end namespace detran


#endif /* WITHINGROUPACCELERATION_HH_ */
