//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EnergyIndependentEigenOperator.cc
 *  @brief EnergyIndependentEigenOperator class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EnergyIndependentEigenOperator.hh"
#include <cstring>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
EnergyIndependentEigenOperator<D>::
EnergyIndependentEigenOperator(SP_mg_solver mg_solver)
  : Base(this)
  , d_mg_solver(mg_solver)
{
  Require(d_mg_solver);

  // Set the fission source
  d_fissionsource = d_mg_solver->fissionsource();

  // Set the operator size
  Assert(d_mg_solver->state());
  set_size(d_mg_solver->state()->moments_size());

  Ensure(d_fissionsource);
}

//----------------------------------------------------------------------------//
template <class D>
void EnergyIndependentEigenOperator<D>::display() const
{
  std::cout << "ENERGY-INDEPENDENT EIGENVALUE OPERATOR" << std::endl;
  std::cout << "              total size: " << d_m << std::endl;
}

//----------------------------------------------------------------------------//
template <class D>
void EnergyIndependentEigenOperator<D>::multiply(const Vector &x,  Vector &y)
{
  // Copy x to the density.
  memcpy(const_cast<double*>(&d_fissionsource->density()[0]),
         &x[0], d_m*sizeof(double));

  // Setup outer iteration.  This precomputes the group sources.
  d_fissionsource->setup_outer();

  // Solve the multigroup equations.
  d_mg_solver->solve();

  // Update the density.
  d_fissionsource->update();

  // Copy density to y.
  memcpy(&y[0], &d_fissionsource->density()[0], d_m*sizeof(double));
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class EnergyIndependentEigenOperator<_1D>;
template class EnergyIndependentEigenOperator<_2D>;
template class EnergyIndependentEigenOperator<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file EnergyIndependentEigenOperator.cc
//----------------------------------------------------------------------------//
