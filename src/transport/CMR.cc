/*
 * CMR.cc
 *
 *  Created on: May 17, 2012
 *      Author: robertsj
 */

#include "CMR.hh"

namespace detran
{

CMR::CMR(SP_mesh mesh, SP_material material, SP_quadrature quadrature)
  : Acceleration(mesh, material, quadrature)
{
  /* ... */
}

} // end namespace detran


