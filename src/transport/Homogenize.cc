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
                       vec_int      coarsegroup,
                       size_t       dc_weight)
{
  Require(state);
  Require(mesh);
  Require(detran_utilities::vec_sum(coarsegroup) == d_number_groups);
  Require(dc_weight < END_DIFF_COEF_WEIGHTING);

  using detran_material::Material;

  // Fine-to-coarse group map
  vec_int fg_to_cg(d_number_groups, 0);
  size_t fg = 0;
  for (size_t cg = 0; cg < coarsegroup.size(); ++cg)
    for (size_t gg = 0; gg < coarsegroup[cg]; ++gg, ++fg)
      fg_to_cg[fg] = cg;

  // Coarse mesh maps
  const vec_int &mesh_map = mesh->mesh_map(key);
  const vec_int &mat_map  = mesh->mesh_map("MATERIAL");

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
    vec2_dbl sigma_s(number_coarse_cells, vec_dbl(d_number_groups, 0.0));
    vec_dbl diff_coef(number_coarse_cells, 0.0);

    // Loop through fine groups in this coarse group
    for (size_t gg = 0; gg < coarsegroup[cg]; ++gg, ++fg)
    {
      const vec_dbl &phi_fg = state->phi(fg);
      const vec_dbl &J_fg   = current(state, fg, dc_weight);

      for (size_t fi = 0; fi < mesh->number_cells(); ++fi)
      {
        size_t ci = mesh_map[fi]; // edit region index
        size_t m  = mat_map[fi];  // material index
        double pv = phi_fg[fi] * mesh->volume(fi);
        double jv = J_fg[fi]   * mesh->volume(fi);
        vol[ci]        += mesh->volume(fi);
        phi_vol[ci]    += pv;
        cur_vol[ci]    += jv;
        sigma_t[ci]    += pv * d_material->sigma_t(m, fg);
        sigma_a[ci]    += pv * d_material->sigma_a(m, fg);
        sigma_f[ci]    += pv * d_material->sigma_f(m, fg);
        nu_sigma_f[ci] += pv * d_material->nu_sigma_f(m, fg);
        chi[ci]        += mesh->volume(fi) * d_material->chi(m, fg);
        for (size_t gp = 0; gp < d_number_groups; ++gp)
          sigma_s[ci][fg_to_cg[gp]] += pv * d_material->sigma_s(m, gp, fg);
        if (dc_weight == PHI_D or dc_weight == CURRENT_D)
          diff_coef[ci] += jv * d_material->diff_coef(m, fg);
      }
    }

    // Set the new coarse mesh values
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
        }
        cmat->set_chi(ci, cg, chi[ci]/vol[ci]);
        for (size_t cgp = 0; cgp < number_coarse_groups; ++cgp)
          cmat->set_sigma_s(ci, cgp, cg, sigma_s[ci][cgp]/phi_vol[ci]);
        if (dc_weight == PHI_D or dc_weight == CURRENT_D)
          cmat->set_diff_coef(ci, cg, diff_coef[ci]/cur_vol[ci]);
        else
          cmat->set_diff_coef(ci, cg, 1.0/(3.0*cmat->sigma_t(ci, cg)));
      }
    }

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
                       std::string  key,
                       size_t       dc_weight)
{
  vec_int coarsegroup(d_number_groups, 1);
  return homogenize(state, mesh, key, coarsegroup, dc_weight);
}

//---------------------------------------------------------------------------//
const Homogenize::vec_dbl&
Homogenize::current(SP_state state, size_t g, size_t dc_weight) const
{
  if (dc_weight == CURRENT_D)
  {
    Insist(state->store_current(),
           "Current-weighting requires the state keeps the current.");
    return state->current(g);
  }
  return state->phi(g);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Homogenize.cc
//---------------------------------------------------------------------------//
