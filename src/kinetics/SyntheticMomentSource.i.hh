//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SyntheticMomentSource.i.hh
 *  @brief  SyntheticMomentSource inline member definitions
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_SYNTHETICMOMENTSOURCE_I_HH_
#define detran_SYNTHETICMOMENTSOURCE_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline double SyntheticMomentSource::source(const size_t cell,
                                      const size_t group)
{
  // Preconditions
  Require(group < d_source.size());
  Require(cell  < d_source[0].size());

  return d_source[group][cell];
}

//---------------------------------------------------------------------------//
inline double SyntheticMomentSource::source(const size_t cell,
                                            const size_t group,
                                            const size_t angle)
{
  // Preconditions
  Require(group < d_source.size());
  Require(cell  < d_source[0].size());

  // Norm is 1/4pi or 1/2
  return d_source[group][cell] * d_norm;
}

//---------------------------------------------------------------------------//
inline void SyntheticMomentSource::build(const double dt,
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
  double a_0 = bdf_coefs[order-1][0];

  // Clear the source
  for (size_t g = 0; g < d_material->number_groups(); ++g)
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
      d_source[g][cell] = 0.0;

  // Add all backward terms
  for (size_t j = 0; j < order; ++j)
  {
    Assert(states[j]);
    if (size_precursor)
    {
      Assert(precursors[j]);
    }

    // Skip the first entry, which is for the (n+1) term
    double a_j = bdf_coefs[order-1][j + 1];

    for (size_t g = 0; g < d_material->number_groups(); ++g)
    {

      // Add the flux term
      double phi_factor = a_j / dt / d_material->velocity(g);
      for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
      {
//        std::cout << " phi factor=" << phi_factor
//                  << " phi = " << states[j]->phi(g)[cell]
//                  << " source = " << d_source[g][cell] << std::endl;
        d_source[g][cell] += phi_factor * states[j]->phi(g)[cell];
      }

      // Add the precursor concentration, if applicable
      if (size_precursor)
      {
        for (size_t i = 0; i < np; ++i)
        {
          double C_factor = a_j * d_material->lambda(i) /
                            (a_0 + dt * d_material->lambda(i));
          for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
          {

            d_source[g][cell] += C_factor *
                                 d_material->chi_d(mt[cell], i, g) *
                                 precursors[j]->C(i)[cell];
//            std::cout << " C=" << precursors[j]->C(i)[cell]
//                      << " chid = " << d_material->chi_d(mt[cell], i, g)
//                      << " q= " << d_source[g][cell]
//                      << std::endl;
          }
        }
      }

    } // end groups
  } // end backward terms
//  std::cout << " source = " << this->source(0, 0) << "  " << this->source(1, 0) << std::endl;
//  std::cout << " source = " << this->source(0, 0, 0) << "  " << this->source(1, 0, 0) << " " << dt << std::endl;
}

} // end namespace detran

#endif // detran_SYNTHETICMOMENTSOURCE_I_HH_

//---------------------------------------------------------------------------//
//              end of file SyntheticMomentSource.i.hh
//---------------------------------------------------------------------------//
