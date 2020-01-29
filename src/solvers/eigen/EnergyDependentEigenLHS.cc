//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EnergyDependentEigenLHS.cc
 *  @brief EnergyDependentEigenLHS class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EnergyDependentEigenLHS.hh"
#include "solvers/mg/MGSolverGMRES.hh"
#include <cstring>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
EnergyDependentEigenLHS<D>::
EnergyDependentEigenLHS(SP_mg_solver mg_solver)
  : Base(this)
  , d_mg_solver(mg_solver)
{
  MGSolverGMRES<D>* mgs =
    dynamic_cast<MGSolverGMRES<D>*>(&(*d_mg_solver->solver()));
  Insist(mgs, "Need MG-GMRES to extract the MG operator.");

  MGTransportOperator<D> &A = *mgs->get_operator();

  // Create scatter-free fission operator
  d_fission =
      new MGScatterFissionOperator(d_mg_solver->input(),
                                   d_mg_solver->material(),
                                   d_mg_solver->mesh(),
                                   A.sweepsource()->get_scatter_source(),
                                   d_mg_solver->fissionsource(),
                                   0, // no cutoff for eigenvalue problems
                                   MGScatterFissionOperator::FISSIONONLY,
                                   d_mg_solver->adjoint());

  // Create the sweep operator get_scatter_source
  d_sweep = new MGSweepOperator<D>(d_mg_solver->state(),
                                   d_mg_solver->boundary(),
                                   A.sweeper(),
                                   A.sweepsource(),
                                   0, // no cutoff for eigenvalue problems
                                   d_mg_solver->adjoint());

  // Operator size
  set_size(d_sweep->number_rows());

  // Size of the moment portion
  d_moments_size = A.moments_size() * d_mg_solver->material()->number_groups();

}

//----------------------------------------------------------------------------//
template <class D>
void EnergyDependentEigenLHS<D>::display() const
{
  std::cout << "ENERGY-INDEPENDENT EIGENVALUE OPERATOR" << std::endl;
  std::cout << "              total size: " << d_m << std::endl;
}

//----------------------------------------------------------------------------//
template <class D>
void EnergyDependentEigenLHS<D>::multiply(const Vector &x,  Vector &y)
{
  // D*inv(L)*M * X * F'

  Vector z(x.size(), 0.0);

  // Shorten z and x and apply fission matrix:   z <-- XF' * x
  Vector X(d_moments_size, const_cast<double*>(&x[0]));
  Vector Z(d_moments_size, &z[0]);
  d_fission->multiply(X, Z);

  // Sweep
  d_sweep->multiply(z, y);
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class EnergyDependentEigenLHS<_1D>;
template class EnergyDependentEigenLHS<_2D>;
template class EnergyDependentEigenLHS<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file EnergyDependentEigenLHS.cc
//----------------------------------------------------------------------------//
