//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenArnoldi.i.hh
 *  @brief EigenArnoldi inline member definitions.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENARNOLDI_I_HH_
#define detran_EIGENARNOLDI_I_HH_

#include "Warning.hh"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include "StdOutUtils.hh"

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
void EigenArnoldi<D>::solve()
{
  std::cout << "Starting EIGEN." << std::endl;

  // Initialize the fission density, copy it to the initial guess,
  // and normalize.
  d_fissionsource->initialize();
  memcpy(&(*d_x0)[0],
         const_cast<double*>(&d_fissionsource->density()[0]),
         d_operator->number_rows()*sizeof(double));
  d_x0->scale(1.0 / d_x0->norm(callow::L2));

  // Solve the problem
  d_eigensolver->solve(d_x, d_x0);

  // Copy the result into the fission source.
  memcpy(const_cast<double*>(&d_fissionsource->density()[0]),
         &(*d_x)[0], d_operator->number_rows()*sizeof(double));

  // To retrieve the correct flux moments, we need to do
  // one more solve with this new density.  This could
  // be a switched feature.
  d_fissionsource->setup_outer();
  d_mg_solver->solve();

  d_state->set_eigenvalue(d_eigensolver->eigenvalue());

  std::cout << "EIGEN done. keff = " << d_eigensolver->eigenvalue()
            <<  std::endl;

}

} // end namespace detran

#endif /* detran_EIGENARNOLDI_I_HH_ */

//----------------------------------------------------------------------------//
//              end of EigenArnoldi.i.hh
//----------------------------------------------------------------------------//

