//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SyntheticDiscreteSource.i.hh
 *  @brief  SyntheticDiscreteSource.i
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_SYNTHETICDISCRETESOURCE_I_HH_
#define detran_SYNTHETICDISCRETESOURCE_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline double SyntheticDiscreteSource::source(const size_t cell,
                                              const size_t group)
{
  // Preconditions
  Require(group < d_source.size());
  Require(cell  < d_source[group][0].size());

  // Integrate over angle
  double value = 0;
  for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
  {
    for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
    {
      size_t angle = d_quadrature->index(o, a);
      value += d_source[group][angle][cell] * d_quadrature->weight(a);
    }
  }
  return value;
}

//---------------------------------------------------------------------------//
inline double SyntheticDiscreteSource::source(const size_t cell,
                                              const size_t group,
                                              const size_t angle)
{
  // Preconditions
  Require(group < d_source.size());
  Require(angle < d_source[group].size());
  Require(cell  < d_source[group][angle].size());

  return d_source[group][angle][cell];
}

//---------------------------------------------------------------------------//
inline void SyntheticDiscreteSource::build(const double dt,
                                           const vec_states &states,
                                           const vec_precursors &precursors,
                                           const size_t order)
{
  // Preconditions
  Require(order > 0);
  Require(order <= 6);
  // Ensure the state size is consistent with the order requested
  size_t size_state = states.size();
  Require(states.size() > 0);
  Require(states.size() >= order);
  // If the precursors are present, they must be the same size as the state
  size_t size_precursor = precursors.size();
  Require((size_precursor == size_state) || (size_precursor == 0));

  // Number of precursors
  size_t np = 0;
  if (size_precursor) np = precursors[0]->number_precursor_groups();

  // Material map
  const detran_utilities::vec_int &mt = d_mesh->mesh_map("MATERIAL");

  // Leading coefficient
  double a_0 = bdf_coefs[order - 1][0];

  // Clear the source
  for (size_t g = 0; g < d_material->number_groups(); ++g)
    for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
      for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
          d_source[g][d_quadrature->index(o, a)][cell] = 0.0;

  // Add all backward terms
  for (size_t j = 0; j < order; ++j)
  {
    Assert(states[j]);
    if (size_precursor)
    {
      Assert(precursors[j]);
    }

    // Skip the first entry, which is for the (n+1) term
    double a_j = bdf_coefs[order - 1][j + 1];

    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {
      double psi_factor = a_j / dt / d_material->velocity(g);

      for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
      {
        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {

          size_t angle = d_quadrature->index(o, a);

          // Add flux
          for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
          {
            d_source[g][angle][cell] +=
              psi_factor * states[j]->psi(g, o, a)[cell];
            //std::cout << " j = " << j << " psi = " << states[j]->psi(g, o, a)[cell] << std::endl;
          }

          // Add the precursor concentration, if applicable
          if (size_precursor)
          {
            for (size_t i = 0; i < np; ++i)
            {
              double C_factor =  a_j * d_norm * d_material->lambda(i) /
                                (a_0 + dt * d_material->lambda(i));
//              std::cout << " Cfactor=" << C_factor
//                        << " den=" << a_0 + dt * d_material->lambda(i)
//                        << " num=" << a_j * d_norm * d_material->lambda(i)
//                        << " aj=" << a_j
//                        << " lam=" << d_material->lambda(i)
//                        << " nrm=" << d_norm << std::endl;
              for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
              {
//                std::cout << " C=" <<precursors[j]->C(i)[cell]
//                          << " chid = " << d_material->chi_d(mt[cell], i, g)
//                          << std::endl;
                d_source[g][angle][cell] += C_factor *
                  d_material->chi_d(mt[cell], i, g) *
                    precursors[j]->C(i)[cell];
              }
            }
          }

        } // end angle
      } // end octant
    } // end groups
  } // end backward terms

}

} // end namespace detran

#endif // detran_SYNTHETICDISCRETESOURCE_I_HH_

//---------------------------------------------------------------------------//
//              end of file SyntheticDiscreteSource.i.hh
//---------------------------------------------------------------------------//
