//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   FissionSource.i.hh
 *  @author robertsj
 *  @date   Oct 14, 2012
 *  @brief  FissionSource inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef detran_FISSIONSOURCE_I_HH_
#define detran_FISSIONSOURCE_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline void FissionSource::setup_outer(const double scale)
{
  d_scale = scale;
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      d_source[g][cell] = d_scale * d_density[cell] *
                          d_material->chi(mat_map[cell], g);
    }
  }
}


//---------------------------------------------------------------------------//
inline void FissionSource::update()
{
  d_density.assign(d_density.size(), 0.0);
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for(size_t g = 0; g < d_number_groups; g++)
  {
    State::moments_type phi = d_state->phi(g);
    for (size_t cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      d_density[cell] += phi[cell] *
                         d_material->nu_sigma_f(mat_map[cell], g);
    }
  }
}

//---------------------------------------------------------------------------//
inline const State::moments_type& FissionSource::source(const size_t g)
{
  Require(g < d_number_groups);
//  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
//  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
//  {
//    d_source[cell] = d_scale * d_density[cell] *
//                     d_material->chi(mat_map[cell], g);
//  }
  return d_source[g];
}

//---------------------------------------------------------------------------//
inline const State::moments_type& FissionSource::density()
{
  return d_density;
}

//---------------------------------------------------------------------------//
inline void FissionSource::
build_within_group_source(const size_t g,
                          const State::moments_type &phi,
                          State::moments_type &source)
{
  Require(g < d_material->number_groups());
  Require(phi.size() == source.size());

  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    source[cell] += phi[cell] * d_scale *
                    d_material->chi(mat_map[cell], g) *
                    d_material->nu_sigma_f(mat_map[cell], g);
  }
}

//---------------------------------------------------------------------------//
inline void FissionSource::
build_in_fission_source(const size_t g,
                        moments_type &source)
{
  Require(g < d_material->number_groups());

  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for (size_t gp = 0; gp < d_material->number_groups(); ++gp)
  {
    if (gp == g) continue;
    const moments_type &phi = d_state->phi(gp);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      source[cell] += phi[cell] * d_scale *
                      d_material->chi(mat_map[cell], g) *
                      d_material->nu_sigma_f(mat_map[cell], gp);
    }
  }

}

//---------------------------------------------------------------------------//
inline void FissionSource::
build_total_group_source(const size_t g,
                         const State::vec_moments_type &phi,
                         State::moments_type &source)
{
  Require(g < d_material->number_groups());

  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for (size_t gp = 0; gp < d_material->number_groups(); ++gp)
  {
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      source[cell] += phi[gp][cell] * d_scale *
                      d_material->chi(mat_map[cell], g) *
                      d_material->nu_sigma_f(mat_map[cell], gp);
    }
  }
}

} // namespace detran

#endif /* detran_FISSIONSOURCE_I_HH_ */
