//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGTransportOperator.cc
 *  @author robertsj
 *  @date   Oct 31, 2012
 *  @brief  MGTransportOperator class definition.
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
  , d_moments_size(0)
  , d_boundary_size(0)
  , d_lower(lower)
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

  // Number of active groups
  d_number_groups = d_state->number_groups() - d_lower;

  // Set the operator size
  set_size(d_number_groups * (d_moments_size + d_boundary_size));

}

//---------------------------------------------------------------------------//
template <class D>
void MGTransportOperator<D>::display() const
{
  std::cout << "WITHIN-GROUP TRANSPORT OPERATOR" << std::endl;
  std::cout << "     total size: " << d_m << std::endl;
  std::cout << "   moments size: " << d_moments_size << std::endl;
  std::cout << "  boundary size: " << d_boundary_size << std::endl;
}

//---------------------------------------------------------------------------//
template <class D>
void MGTransportOperator<D>::multiply(const Vector &x,  Vector &y)
{
  using std::cout;
  using std::endl;

  State::vec_moments_type phi_original(d_state->all_phi());
  State::vec_moments_type phi_update(d_state->all_phi());

  // Fill temporary multigroup flux vector
  for (int g = d_lower; g < d_number_groups; g++)
  {
    int offset = (g - d_lower) * (d_moments_size + d_boundary_size);
    for (int i = 0; i < d_moments_size; i++)
    {
      phi_original[g][i] = x[i + offset];
      phi_update[g][i]   = x[i + offset];
    }
  }

  // \todo Switch indexing so that moments contiguous and
  //       boundaries are tacked on the end.  This will make
  //       diffusion preconditioning easier (I think).
  //

  // Sweep each applicable group
  for (int g = d_lower; g < d_number_groups; g++)
  {

    int g_index  = g - d_lower;
    int g_size   = d_moments_size + d_boundary_size;
    // moment offset
    int m_offset = g_index * g_size;
    // boundary offset
    int b_offset = m_offset + d_moments_size;

    // Reset the source and place the original outgoing boundary flux.
    d_boundary->clear(g);

    if (d_boundary->has_reflective())
    {
      // Set the incident boundary flux.
      d_boundary->psi(g, const_cast<double*>(&x[0]) + b_offset,
                      BoundaryBase<D>::IN, BoundaryBase<D>::SET, true);
    }

    // Reset the source to zero.
    d_sweepsource->reset();
    d_sweepsource->build_total_scatter(g, d_lower, phi_original);

    // Set the sweeper and sweep.
    d_sweeper->setup_group(g);
    d_sweeper->sweep(phi_update[g]);

    // Update outgoing vector
    {
      // Assign the moment values.
      for (int i = 0; i < d_moments_size; i++)
      {
        y[i + m_offset] = phi_original[g][i] - phi_update[g][i];
      }
      // Assign boundary fluxes, if applicable
      if (d_boundary->has_reflective())
      {
        // Update the boundary and fetch.
        d_boundary->update(g);
        State::angular_flux_type psi_update(d_boundary_size, 0.0);
        d_boundary->psi(g, &psi_update[0],
                        BoundaryBase<D>::IN, BoundaryBase<D>::GET, true);

        // Add the boundary values.
        for (int a = 0; a < d_boundary_size; a++)
        {
          y[a + b_offset] = x[a + b_offset] - psi_update[a];
        }
      }
    } // end update

   }

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



