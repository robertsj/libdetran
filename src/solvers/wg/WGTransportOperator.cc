//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  WGTransportOperator.cc
 *  @brief WGTransportOperator member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "WGTransportOperator.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
WGTransportOperator<D>::WGTransportOperator(SP_state        state,
                                            SP_boundary     boundary,
                                            SP_sweeper      sweeper,
                                            SP_sweepsource  source)
  : Base(this)
  , d_state(state)
  , d_boundary(boundary)
  , d_sweeper(sweeper)
  , d_sweepsource(source)
  , d_moments_size(0)
  , d_boundary_size(0)
  , d_g(0)
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
void WGTransportOperator<D>::set_group(const size_t g)
{
  d_g = g;
}


//---------------------------------------------------------------------------//
template <class D>
void WGTransportOperator<D>::display() const
{
  std::cout << "WITHIN-GROUP TRANSPORT OPERATOR" << std::endl;
  std::cout << "     total size: " << d_m << std::endl;
  std::cout << "   moments size: " << d_moments_size << std::endl;
  std::cout << "  boundary size: " << d_boundary_size << std::endl;
}

//---------------------------------------------------------------------------//
template <class D>
void WGTransportOperator<D>::multiply(const Vector &x,  Vector &y)
{

  // Fill a temporary flux vector with the Krylov vector.
  State::moments_type phi(d_moments_size, 0.0);
  for (int i = 0; i < d_moments_size; i++) phi[i] = x[i];

  // Reset the boundary, and insert the boundary portion of the Krylov
  // vector, if applicable.
  d_boundary->clear(d_g);
  if (d_boundary->has_reflective())
    d_boundary->psi(d_g, const_cast<double*>(&x[0]) + d_moments_size,
                    BoundaryBase<D>::IN, BoundaryBase<D>::SET, true);

  // Reset the source and build the scattering source (with no
  // fixed source contributions)
  d_sweepsource->reset();
  d_sweepsource->build_within_group_scatter(d_g, phi);

  // Sweep.  This gives X <-- D*inv(L)*M*S*X
  d_sweeper->sweep(phi);

  // Assign the moment values.  This gives X <- (I-D*inv(L)*M*S)*X
  for (int i = 0; i < d_moments_size; i++)
    y[i] = x[i] - phi[i];

  if (d_boundary->has_reflective())
  {
    // Update the boundary (redirect outgoing as incident)
    d_boundary->update(d_g);

    // Extract the incident boundary
    State::moments_type psi_update(d_boundary_size, 0.0);
    d_boundary->psi(d_g, &psi_update[0],
                    BoundaryBase<D>::IN, BoundaryBase<D>::GET, true);

    // Add the boundary values
    for (int a = 0; a < d_boundary_size; ++a)
      y[a + d_moments_size] = x[a + d_moments_size] - psi_update[a];
  }

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class WGTransportOperator<_1D>;
template class WGTransportOperator<_2D>;
template class WGTransportOperator<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file WGTransportOperator.cc
//---------------------------------------------------------------------------//
