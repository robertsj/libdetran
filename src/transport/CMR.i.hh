/*
 * CMR.i.hh
 *
 *  Created on: May 18, 2012
 *      Author: robertsj
 */

#ifndef CMR_I_HH_
#define CMR_I_HH_

namespace detran
{

// 3D
template <class D>
inline void CMR<D>::tally(int i, int j, int k, int o, int a, face_flux_type psi)
{
  return;
}

template <>
inline void CMR<_1D>::tally(int i, int j, int k, int o, int a, face_flux_type psi)
{
  Require(j == 0);
  Require(k == 0);
  Require(o == 0 || o == 1);

  double mu = b_quadrature->mu(0, a);
  double wt = b_quadrature->weight(a);
  double val = psi * mu * wt;

  if (o == 0)
  {
    d_J_pos[0][i+1] += val;
  }
  else
  {
    d_J_neg[0][i-1] += val;
  }

  return;
}
} // end namespace detran


#endif /* CMR_I_HH_ */
