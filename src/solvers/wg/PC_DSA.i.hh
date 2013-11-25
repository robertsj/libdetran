//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  PC_DSA.i.hh
 *  @brief PC_DSA inline member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_PC_DSA_I_HH_
#define detran_PC_DSA_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline void PC_DSA::apply(Vector &V_in, Vector &V_out)
{
  // Currently, DSA is only used on the flux moments,
  // not on the boundaries.
  size_t size_moments = d_mesh->number_cells();

  // Temporary vectors
  State::moments_type S_times_V(size_moments, 0.0);
  State::moments_type V(size_moments, 0.0);
  Vector tmp(size_moments, 0.0);
  Vector V_out_trunc(size_moments, 0.0);

  V_out.set(0.0);

  //-------------------------------------------------------------------------//
  // BEGIN WITH V0.  (Given x, copy into xv for manipulation)
  //-------------------------------------------------------------------------//

  for (int i = 0; i < size_moments; i++)
    V[i] = V_in[i];

  //-------------------------------------------------------------------------//
  // OPERATE: V1 <-- S*V0.
  //-------------------------------------------------------------------------//

  d_scattersource->build_within_group_source(d_group, V, S_times_V);

  // Copy to vector
  for (int i = 0; i < size_moments; i++)
    tmp[i] = S_times_V[i];

  //-------------------------------------------------------------------------//
  // OPERATE: V2 <-- inv(C)*V1 = inv(C)*S*V0.
  //-------------------------------------------------------------------------//

  // V_out <-- C\(S*V_in)
//  if (d_solver[d_group]->preconditioner())
//    std::cout << d_solver[d_group]->preconditioner()->name() << std::endl;
  d_solver[d_group]->solve(tmp, V_out_trunc);

  //-------------------------------------------------------------------------//
  // OPERATE: V3 <-- V0 + V2 = V0 + inv(C)*S*V0 = (I + inv(C)*S)*V0
  //-------------------------------------------------------------------------//

  V_out.copy(V_in);

  for (int i = 0; i < size_moments; ++i)
    V_out[i] += V_out_trunc[i];

}

} // end namespace detran

#endif // detran_PC_DSA_I_HH_

//---------------------------------------------------------------------------//
//              end of file PreconditionerWG.i.hh
//---------------------------------------------------------------------------//
