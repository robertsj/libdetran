//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Homogenize.cc
 *  @brief  Homogenize 
 *  @author Jeremy Roberts
 *  @date   Mar 26, 2013
 */
//---------------------------------------------------------------------------//

#include "Homogenize.hh"
#include "utilities/MathUtilities.hh"
#include "material/Material.hh"

namespace detran
{

//---------------------------------------------------------------------------//
Homogenize::Homogenize(SP_material material)
  : d_material(material)
{
  d_number_groups = d_material->number_groups();
}

//---------------------------------------------------------------------------//
Homogenize::SP_material
Homogenize::homogenize(SP_state     state,
                      SP_mesh      mesh,
                      std::string  key,
                      vec_int      coarsegroup)
{
  Require(state);
  Require(mesh);
  Require(detran_utilities::vec_sum(coarsegroup) == d_number_groups);

  using detran_material::Material;

  // Fine-to-coarse group map
  vec_int fg_to_cg(d_number_groups, 0);
  size_t fg = 0;
  for (size_t cg = 0; cg < coarsegroup.size(); ++cg)
    for (size_t gg = 0; gg < coarsegroup[cg]; ++gg, ++fg)
      fg_to_cg[fg] = cg;

  // Coarse mesh map
  const vec_int &mesh_map = mesh->mesh_map(key);

  // Maximum number of unique coarse cells, i.e. materials
  size_t number_coarse_cells = detran_utilities::vec_max(mesh_map) + 1;

  // Number of coarse groups
  size_t number_coarse_groups = coarsegroup.size();

  // Create new material
  SP_material cmat =
    Material::Create(number_coarse_cells, number_coarse_groups);

  // Condense
  fg = 0;
  for (size_t cg = 0; cg < number_coarse_groups; ++cg)
  {
    // Zero out temporaries
    vec_dbl vol(number_coarse_cells, 0.0);
    vec_dbl phi_vol(number_coarse_cells, 0.0);
    vec_dbl cur_vol(number_coarse_cells, 0.0);
    vec_dbl sigma_t(number_coarse_cells, 0.0);
    vec_dbl sigma_a(number_coarse_cells, 0.0);
    vec_dbl nu_sigma_f(number_coarse_cells, 0.0);
    vec_dbl sigma_f(number_coarse_cells, 0.0);
    vec_dbl nu(number_coarse_cells, 0.0);
    vec_dbl chi(number_coarse_cells, 0.0);
    vec2_dbl sigma_s(number_coarse_cells,
                     vec_dbl(d_number_groups, 0.0));
    vec_dbl diff_coef(number_coarse_cells, 0.0);

    // Loop through fine groups in this coarse group
    for (size_t gg = 0; gg < coarsegroup[cg]; ++gg, ++fg)
    {
      vec_dbl &phi_fg = state->phi(fg);
      for (size_t fi = 0; fi < mesh->number_cells(); ++fi)
      {
        size_t ci = mesh_map[fi];
        double pv = phi_fg[fi] * mesh->volume(fi);
        vol[ci]        += mesh->volume(fi);
        phi_vol[ci]    += pv;
        sigma_t[ci]    += pv * d_material->sigma_t(ci, fg);
        sigma_a[ci]    += pv * d_material->sigma_a(ci, fg);
        sigma_f[ci]    += pv * d_material->sigma_f(ci, fg);
        nu_sigma_f[ci] += pv * d_material->nu_sigma_f(ci, fg);
        chi[ci]        += mesh->volume(fi) * d_material->chi(ci, fg);
        for (size_t gp = 0; gp < d_number_groups; ++gp)
          sigma_s[ci][fg_to_cg[gp]] += pv * d_material->sigma_s(ci, gp, fg);
        // \todo Add options for homogenization.  Ideally, current-weighted.
        diff_coef[ci] += pv * d_material->diff_coef(ci, fg);
      }
    } // end energy

    for (size_t ci = 0; ci < number_coarse_cells; ++ci)
    {
      // we might have indices that don't show up in the problem.
      if (sigma_t[ci] > 0.0)
      {
        cmat->set_sigma_t(ci, cg, sigma_t[ci]/phi_vol[ci]);
        cmat->set_sigma_a(ci, cg, sigma_a[ci]/phi_vol[ci]);
        if (sigma_f[ci] > 0.0)
        {
          cmat->set_sigma_f(ci, cg, sigma_f[ci]/phi_vol[ci]);
          cmat->set_nu_sigma_f(ci, cg, nu_sigma_f[ci]/phi_vol[ci]);
          cmat->set_nu(ci, cg, nu_sigma_f[ci]/sigma_f[ci]);
          cmat->set_chi(ci, cg, chi[ci]/vol[ci]);
        }

        for (size_t cgp = 0; cgp < number_coarse_groups; ++cgp)
          cmat->set_sigma_s(ci, cgp, cg, sigma_s[ci][cgp]/phi_vol[ci]);
        cmat->set_diff_coef(ci, cg, diff_coef[ci]/phi_vol[ci]);
      }
    } // end space

  } // end coarse energy

  // normalize chi
  for (int m = 0; m < number_coarse_cells; ++m)
  {
    double chi_sum = 0.0;
    for (int g = 0; g < number_coarse_groups; ++g)
      chi_sum += cmat->chi(m, g);
    if (chi_sum > 0.0)
      for (int g = 0; g < number_coarse_groups; ++g)
        cmat->set_chi(m, g, cmat->chi(m, g)/chi_sum);
  }

  cmat->finalize();
  return cmat;
}

//---------------------------------------------------------------------------//
Homogenize::SP_material
Homogenize::homogenize(SP_state     state,
                       SP_mesh      mesh,
                       std::string  key)
{
  vec_int coarsegroup(d_number_groups, 1);
  return homogenize(state, mesh, key, coarsegroup);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Homogenize.cc
//---------------------------------------------------------------------------//
