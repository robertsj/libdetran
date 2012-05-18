/*
 * CMR.hh
 *
 *  Created on: May 17, 2012
 *      Author: robertsj
 */

#ifndef CMR_HH_
#define CMR_HH_

// Detran
#include "Acceleration.hh"

// System

namespace detran
{



/*!
 *  \class CMR
 *  \brief Coarse Mesh Rebalance
 *
 *  Note, this applies only to the within group problem.
 */
class CMR : public Acceleration
{

public:

  /// \name Useful Typedefs
  // \{

  typedef SP<CMR>                   SP_acceleration;
  typedef Acceleration              Base;
  typedef Base::SP_acceleration     SP_base;
  typedef Mesh::SP_mesh             SP_mesh;
  typedef Material::SP_material     SP_material;
  typedef Quadrature::SP_quadrature SP_quadrature;
  typedef State::SP_state           SP_state;

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

  /// See \ref Acceleration::tally
  void tally(int i, int j, int k, int o, int a, int g, double psi)
  {
    return;
  }

private:

};


}

#endif /* CMR_HH_ */
