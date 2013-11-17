//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Homogenize.cc
 *  @brief Homogenize member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Homogenize.hh"
#include "material/Material.hh"
#include "utilities/MathUtilities.hh"

namespace detran
{

using detran_utilities::range;
using detran_utilities::vec_sum;

//----------------------------------------------------------------------------//
Homogenize::Homogenize(SP_material material, const size_t dc_weight)
  : d_material(material)
  , d_option_dc(dc_weight)
  , d_option_spectrum(END_OPTIONS_SPECTRUM)
{
  Require(dc_weight < END_OPTIONS_DC);
  d_number_groups = d_material->number_groups();
}

//----------------------------------------------------------------------------//
Homogenize::SP_material
Homogenize::homogenize(SP_state            state,
                       SP_mesh             mesh,
                       const std::string  &regionkey,
                       vec_int             coarsegroup,
                       vec_size_t          groups)
{
  Require(state);
  Require(mesh);
  d_state = state;
  d_option_spectrum = FINE_MESH_SPECTRUM;
  // Only if the state stores the current can we current-weight for D.
  if (d_option_dc == CURRENT_D and !d_state->store_current())
    d_option_dc = PHI_D;
  return homogenize(mesh, regionkey, coarsegroup, groups);
}

//----------------------------------------------------------------------------//
Homogenize::SP_material
Homogenize::homogenize(const vec2_dbl    &spectrum,
                       const std::string &spectrumkey,
                       SP_mesh            mesh,
                       const std::string &regionkey,
                       vec_int            coarsegroup,
                       vec_size_t         groups)
{
  Require(spectrum.size() == d_material->number_groups());
  Require(mesh);
  d_spectrum = spectrum;
  d_spectrum_map = mesh->mesh_map(spectrumkey);
  d_option_spectrum = REGION_SPECTRUM;
  return homogenize(mesh, regionkey, coarsegroup, groups);
}

//----------------------------------------------------------------------------//
Homogenize::SP_material
Homogenize::homogenize(SP_mesh            mesh,
                       const std::string &regionkey,
                       vec_int            coarsegroup,
                       vec_size_t         groups)
{
  if (!groups.size()) groups = range<size_t>(0, d_number_groups);
  if (!coarsegroup.size()) coarsegroup = vec_int(groups.size(), 1);
  Require(groups.size() > 0);
  Require(vec_sum(coarsegroup) == groups.size());

  // Reverse groups if adjoint
  if (*(groups.begin()) > *(groups.end()-1))
    std::reverse(groups.begin(), groups.end());

  // Fine-to-coarse group map
  vec_int fg_to_cg(d_number_groups, 0);
  size_t fg = 0;
  for (size_t cg = 0; cg < coarsegroup.size(); ++cg)
    for (int gg = 0; gg < coarsegroup[cg]; ++gg, ++fg)
      fg_to_cg[fg] = cg;

  // Region mesh maps
  const vec_int &reg_map = mesh->mesh_map(regionkey);
  const vec_int &mat_map = mesh->mesh_map("MATERIAL");

  // Maximum number of unique coarse cells.  Note, if that map has missing
  // indices, then empty materials will be produced.  If we changed the
  // material class to hold materials separately, we could just have empty
  // pointers where a material doesn't exist.  In practice, it shouldn't
  // matter (outside of storing extra zeros).
  size_t num_coarse_cells = detran_utilities::vec_max(reg_map) + 1;

  // Number of coarse groups
  size_t num_coarse_groups = coarsegroup.size();

  // Create new material
  SP_material cmat =
      detran_material::Material::Create(num_coarse_cells, num_coarse_groups);

  // Size the averaged flux
  d_coarse_mesh_flux.resize(num_coarse_groups, vec_dbl(num_coarse_cells, 0.0));

  // Loop over all coarse groups
  fg = 0;
  for (size_t cg = 0; cg < num_coarse_groups; ++cg)
  {
    // Zero out temporaries
    vec_dbl  vol(num_coarse_cells, 0.0);
    vec_dbl  phi_vol(num_coarse_cells, 0.0);
    vec_dbl  cur_vol(num_coarse_cells, 0.0);
    vec_dbl  sigma_t(num_coarse_cells, 0.0);
    vec_dbl  sigma_a(num_coarse_cells, 0.0);
    vec_dbl  nu_sigma_f(num_coarse_cells, 0.0);
    vec_dbl  sigma_f(num_coarse_cells, 0.0);
    vec_dbl  nu(num_coarse_cells, 0.0);
    vec_dbl  chi(num_coarse_cells, 0.0);
    vec_dbl  fd(num_coarse_cells, 0.0);
    vec2_dbl sigma_s(num_coarse_cells, vec_dbl(d_number_groups, 0.0));
    vec_dbl  diff_coef(num_coarse_cells, 0.0);

    // Loop through fine groups in this coarse group
    for (int gg = 0; gg < coarsegroup[cg]; ++gg, ++fg)
    {
      // Actual fine group
      size_t group = fg + groups[0];

      // Loop over cells
      for (size_t fi = 0; fi < mesh->number_cells(); ++fi)
      {
        size_t ci = reg_map[fi]; // edit region index
        size_t m  = mat_map[fi]; // material index
        double pv = spectrum(fg, fi) * mesh->volume(fi);
        double jv = current(fg, fi)  * mesh->volume(fi);
        vol[ci]        += mesh->volume(fi);
        phi_vol[ci]    += pv;
        cur_vol[ci]    += jv;
        sigma_t[ci]    += pv * d_material->sigma_t(m, group);
        sigma_a[ci]    += pv * d_material->sigma_a(m, group);
        sigma_f[ci]    += pv * d_material->sigma_f(m, group);
        nu_sigma_f[ci] += pv * d_material->nu_sigma_f(m, group);
        // Volume weighting chi with subsequent normalization to unity
        //chi[ci]        += mesh->volume(fi) * d_material->chi(m, group);
        double temp_fd = 0.0;
        for (size_t gp = 0; gp < d_number_groups; ++gp)
          temp_fd = spectrum(gp, fi) * d_material->nu_sigma_f(m, gp) * mesh->volume(fi);
        if (gg==0) fd[ci] += temp_fd;
        chi[ci] += temp_fd * d_material->chi(m, group);
        for (size_t gp = 0; gp < d_number_groups; ++gp)
          sigma_s[ci][fg_to_cg[gp]] += pv * d_material->sigma_s(m, gp, group);
        if (d_option_dc == PHI_D || d_option_dc == CURRENT_D)
          diff_coef[ci] += jv * d_material->diff_coef(m, group);
      }
    }

    // Set the new coarse mesh values
    for (size_t ci = 0; ci < num_coarse_cells; ++ci)
    {
      // We might have indices that don't show up in the problem.
      if (sigma_t[ci] > 0.0)
      {
        cmat->set_sigma_t(ci, cg, sigma_t[ci]/phi_vol[ci]);
        cmat->set_sigma_a(ci, cg, sigma_a[ci]/phi_vol[ci]);
        if (sigma_f[ci] > 0.0)
        {
          cmat->set_sigma_f(ci, cg, sigma_f[ci]/phi_vol[ci]);
          cmat->set_nu_sigma_f(ci, cg, nu_sigma_f[ci]/phi_vol[ci]);
          cmat->set_nu(ci, cg, nu_sigma_f[ci]/sigma_f[ci]);
          cmat->set_chi(ci, cg, chi[ci] / fd[ci]);
        }
        //cmat->set_chi(ci, cg, chi[ci]/vol[ci]);
        for (size_t cgp = 0; cgp < num_coarse_groups; ++cgp)
          cmat->set_sigma_s(ci, cgp, cg, sigma_s[ci][cgp]/phi_vol[ci]);
        if (d_option_dc == PHI_D || d_option_dc == CURRENT_D)
          cmat->set_diff_coef(ci, cg, diff_coef[ci]/cur_vol[ci]);
        else
          cmat->set_diff_coef(ci, cg, 1.0/(3.0*cmat->sigma_t(ci, cg)));
      }

      d_coarse_mesh_flux[cg][ci] = phi_vol[ci] / vol[ci];
    }

  } // end coarse energy

  // normalize chi
//  for (size_t m = 0; m < num_coarse_cells; ++m)
//  {
//    double chi_sum = 0.0;
//    for (size_t g = 0; g < num_coarse_groups; ++g)
//      chi_sum += cmat->chi(m, g);
//    if (chi_sum > 0.0)
//      for (size_t g = 0; g < num_coarse_groups; ++g)
//        cmat->set_chi(m, g, cmat->chi(m, g)/chi_sum);
//  }

  cmat->finalize();
  return cmat;
}

//----------------------------------------------------------------------------//
void Homogenize::set_option_dc(const size_t option_dc)
{
  Require(option_dc < END_OPTIONS_DC);
  d_option_dc = option_dc;
}

//----------------------------------------------------------------------------//
double Homogenize::current(const size_t g, const size_t cell) const
{
  Require(d_option_spectrum < END_OPTIONS_SPECTRUM);
  Require(g < d_material->number_groups());
  double value = 0.0;
  if (d_option_spectrum == FINE_MESH_SPECTRUM)
  {
    if (d_option_dc == CURRENT_D)
      value = d_state->current(g)[cell];
    else
      value = d_state->phi(g)[cell];
  }
  else if (d_option_spectrum == REGION_SPECTRUM)
  {
    value = d_spectrum[g][d_spectrum_map[cell]];
  }
  return value;
}

//----------------------------------------------------------------------------//
double Homogenize::spectrum(const size_t g, const size_t cell) const
{
  Require(d_option_spectrum < END_OPTIONS_SPECTRUM);
  Require(g < d_material->number_groups());
  double value = 0.0;
  if (d_option_spectrum == FINE_MESH_SPECTRUM)
  {
    value = d_state->phi(g)[cell];
  }
  else if (d_option_spectrum == REGION_SPECTRUM)
  {
    value = d_spectrum[g][d_spectrum_map[cell]];
  }
  return value;
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file Homogenize.cc
//----------------------------------------------------------------------------//
