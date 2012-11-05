//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGTransportOperator.cc
 *  @brief  MGTransportOperator member definitions
 *  @author robertsj
 *  @date   Oct 31, 2012
 */
//---------------------------------------------------------------------------//

#include "MGTransportOperator.hh"

namespace detran
{

template <class D>
MGTransportOperator<D>::MGTransportOperator(SP_state        state,
                                            SP_boundary     boundary,
                                            SP_sweeper      sweeper,
                                            SP_sweepsource  source,
                                            size_t          lower)
  : Base(this)
  , d_state(state)
  , d_boundary(boundary)
  , d_sweeper(sweeper)
  , d_sweepsource(source)
  , d_lower(lower)
  , d_moments_size(0)
  , d_boundary_size(0)
{
  // Preconditions
  Require(d_state);
  Require(d_boundary);
  Require(d_sweeper);
  Require(d_sweepsource);
  Require(d_lower < d_state->number_groups());

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

  // Number of active groups
  d_number_groups = d_state->number_groups() - d_lower;

  // Set the operator size
  set_size(d_number_groups * (d_moments_size + d_boundary_size));

}

//---------------------------------------------------------------------------//
template <class D>
void MGTransportOperator<D>::display() const
{
  std::cout << "MULTI-GROUP TRANSPORT OPERATOR" << std::endl;
  std::cout << "       total size: " << d_m << std::endl;
  std::cout << " number of groups: " << d_number_groups << std::endl;
  std::cout << "     moments size: " << d_moments_size << std::endl;
  std::cout << "    boundary size: " << d_boundary_size << std::endl;
}

//---------------------------------------------------------------------------//
template <class D>
void MGTransportOperator<D>::multiply(const Vector &x,  Vector &y)
{
  using std::cout;
  using std::endl;

  // Fill a temporary flux vector with the Krylov vector.  The Krylov
  // vector is  ordered as follows:
  //   [phi(lower) phi(lower+1) ... phi(G-1) psi(lower+1) ...]
  // where psi is the boundary angular flux, present only if there
  // are reflective conditions
  State::vec_moments_type phi(d_state->all_phi());
  for (int g = d_lower; g < d_number_groups; g++)
  {
    int offset = (g - d_lower) * d_moments_size;
    for (int i = 0; i < d_moments_size; i++)
      phi[g][i] = x[i + offset];
  }

  // sweep each applicable group
  for (int g = d_lower; g < d_number_groups; g++)
  {
    // group index in applicable set
    int g_index  = g - d_lower;
    // moment offset, the starting moment index within the Krylov vector
    int m_offset = g_index * d_moments_size;
    // boundary offset, the starting boundary index within the Krylov vector
    int b_offset = d_number_groups * d_moments_size + g_index * d_boundary_size;

    // reset the source and place the original outgoing boundary flux.
    d_boundary->clear(g);

    if (d_boundary->has_reflective())
    {
      // set the incident boundary flux.
      d_boundary->psi(g, const_cast<double*>(&x[0]) + b_offset,
                      BoundaryBase<D>::IN, BoundaryBase<D>::SET, true);
    }

    // reset the source to zero.
    d_sweepsource->reset();
    d_sweepsource->build_total_scatter(g, d_lower, phi);

    // set the sweeper and sweep.
    d_sweeper->setup_group(g);
    d_sweeper->sweep(phi[g]);

    // update boundary component
    {
      // assign the moment values.
      for (int i = 0; i < d_moments_size; i++)
        y[i + m_offset] = x[i + m_offset] - phi[g][i];

      // assign boundary fluxes, if applicable
      if (d_boundary->has_reflective())
      {
        // update the boundary (redirect outgoing as incident)
        d_boundary->update(g);

        // extract the incident boundary
        State::angular_flux_type psi_update(d_boundary_size, 0.0);
        d_boundary->psi(g, &psi_update[0],
                        BoundaryBase<D>::IN, BoundaryBase<D>::GET, true);

        // add the boundary values.
        for (int a = 0; a < d_boundary_size; a++)
          y[a + b_offset] = x[a + b_offset] - psi_update[a];
      }
    } // end boundary

  } // end groups

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class MGTransportOperator<_1D>;
template class MGTransportOperator<_2D>;
template class MGTransportOperator<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file MGTransportOperator.cc
//---------------------------------------------------------------------------//



