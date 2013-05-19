//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  SweepOperator.cc
 *  @brief SweepOperator
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "SweepOperator.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
SweepOperator<D>::SweepOperator(SP_state        state,
                                SP_boundary     boundary,
                                SP_sweeper      sweeper,
                                SP_sweepsource  source)
  : Base(this)
  , d_state(state)
  , d_boundary(boundary)
  , d_sweeper(sweeper)
  , d_sweepsource(source)
{
  // Preconditions
  Require(d_state);
  Require(d_boundary);
  Require(d_sweeper);
  Require(d_sweepsource);
  Require(d_state->get_mesh());

  // Determine the sizes of the moments.
  d_moments_size = d_state->moments_size();

  // Determine the sizes of any reflected boundary fluxes.
  for (int side = 0; side < 2*D::dimension; side++)
  {
    if (boundary->is_reflective(side))
    {
      // We only need to store only half of the unknowns.
      d_boundary_size += boundary->boundary_flux_size(side)/2;
    }
  }

  // Set the operator size
  set_size(d_moments_size + d_boundary_size);

}

//---------------------------------------------------------------------------//
template <class D>
void SweepOperator<D>::display() const
{
  std::cout << "SWEEP OPERATOR" << std::endl;
  std::cout << "     total size: " << d_m << std::endl;
}

//---------------------------------------------------------------------------//
template <class D>
void SweepOperator<D>::multiply(const Vector &x,  Vector &y)
{

  // Moment size
  size_t size = x.size();

  // Fill a temporary flux vector with the Krylov vector.
  State::moments_type phi(x.size(), 0.0);
  for (int i = 0; i < size; i++) phi[i] = x[i];

  // Sweep.  This gives X <-- D*inv(L)*X
  d_sweeper->sweep(phi);

  // Fill the outgoing vector.  It's time to rethink
  // the state vector.
  for (int i = 0; i < size; i++) y[i] = phi[i];

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class SweepOperator<_1D>;
template class SweepOperator<_2D>;
template class SweepOperator<_3D>;

} // end namespace detran
//---------------------------------------------------------------------------//
//              end of file SweepOperator.cc
//---------------------------------------------------------------------------//
