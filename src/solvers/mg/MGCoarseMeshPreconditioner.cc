//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCoarseMeshPreconditioner.cc
 *  @brief MGCoarseMeshPreconditioner member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGCoarseMeshPreconditioner.hh"
#include "SpectrumGS.hh"
#include "SpectrumFS.hh"
#include "SpectrumPinCell.hh"
#include "transport/Homogenize.hh"
#include "utilities/MathUtilities.hh"

#define COUT(c) std::cout << c << std::endl;

namespace detran
{

//----------------------------------------------------------------------------//
MGCoarseMeshPreconditioner::
MGCoarseMeshPreconditioner(SP_input         input,
                           SP_material      material,
                           SP_mesh          mesh,
                           SP_scattersource ssource,
                           SP_fissionsource fsource,
                           size_t           cutoff,
                           bool             include_fission,
                           bool             adjoint,
                           std::string      name)
  : Base(input, material, mesh, ssource, fsource, cutoff,
         include_fission, adjoint, name)
  , d_condensation_option(CONDENSE_WITH_FLAT_SPECTRUM)
  , d_fine_per_coarse(d_number_active_groups, 1)
  , d_f2c_group_map(d_number_active_groups, 0)
{
  using detran_utilities::vec_sum;

  // In the constructor, things are set that won't change.  Anything
  // that is spectrum-dependent is kept in the separate build routine so
  // that keff-dependent or time-dependent materials can be included in
  // any client-initiated update.

  if (d_input->check("mgpc_condensation_option"))
  {
    d_condensation_option = d_input->get<int>("mgpc_condensation_option");
    Insist(d_condensation_option < END_CONDENSATION_OPTIONS,
           "Invalid MGPC condensation option selected.");
  }

  // Define the active groups of the operator
  int upper = d_adjoint ? -1 : d_number_groups;
  d_groups = detran_utilities::range<size_t>(d_group_cutoff, upper);

  // Build the coarse mesh
  size_t level = 2;
  if (d_input->check("mgpc_coarse_mesh_level"))
  {
    level = d_input->get<int>("mgpc_coarse_mesh_level");
    Assert(level > 0);
  }
  CoarseMesh C(d_mesh, level);
  d_coarsemesh = C.get_coarse_mesh();

  // Define fine to coarse map
  if (d_input->check("mgpc_fine_per_coarse_group"))
  {
    d_fine_per_coarse = d_input->get<vec_int>("mgpc_fine_per_coarse_group");
    Require(vec_sum(d_fine_per_coarse) == d_number_active_groups);
  }

  size_t g = 0;
  for (size_t cg = 0; cg < d_fine_per_coarse.size(); ++cg)
    for (size_t i = 0; i < d_fine_per_coarse[cg]; ++i, ++g)
      d_f2c_group_map[g] = cg;

  d_size_fine   = d_number_active_groups * d_mesh->number_cells();
  d_size_coarse = d_fine_per_coarse.size() * d_coarsemesh->number_cells();
  d_size        = d_size_fine;
  build_restrict();
}

//----------------------------------------------------------------------------//
MGCoarseMeshPreconditioner::~MGCoarseMeshPreconditioner()
{
  /* ... */
}

//----------------------------------------------------------------------------//
void MGCoarseMeshPreconditioner::build(const double keff, SP_state state)
{
  Require(keff != 0.0);

  using detran_utilities::vec_max;
  using detran_utilities::vec_unique;

  try
  {

  //--------------------------------------------------------------------------//
  // SPECTRUM FOR WEIGHTING
  //--------------------------------------------------------------------------//

  if (d_condensation_option == CONDENSE_WITH_STATE)
  {
    THROW("NOT IMPLEMENTED");
    Insist(state, "To build a MGPC using a state requires state != NULL");
    d_state = state;
  }
  else
  {
    // Note, the spectrum must have all groups, even though the material might
    // not be homogenized over all groups. In the future, implement the other
    // ideas, which are in general dependent on time and/or keff.

    size_t number_regions = 0;
    if (d_condensation_option == CONDENSE_WITH_FLAT_SPECTRUM)
    {
      d_spectrum.resize(d_number_groups,
                        vec_dbl(d_material->number_materials(), 1.0));
      d_key = "MATERIAL";
      number_regions = d_material->number_materials();
    }
    else if (d_condensation_option == CONDENSE_WITH_USER_SPECTRUM)
    {
      Insist(d_input->check("mgpc_spectrum"),
             "A user-defined spectrum must be set via mgpc_spectrum");
      Insist(d_input->check("mgpc_spectrum_key"),
             "A user-defined spectrum key must be set via mgpc_spectrum_key");

      vec_dbl spectrum_1d = d_input->get<vec_dbl>("mgpc_spectrum");
      Assert(spectrum_1d.size() % d_number_groups == 0);

      d_key = d_input->get<std::string>("mgpc_spectrum_key");
      Assert(d_mesh->mesh_map_exists(d_key));

      number_regions = spectrum_1d.size() / d_number_groups;
      Assert(number_regions > 0);
      Assert(number_regions == vec_max(d_mesh->mesh_map(d_key))+1);
      size_t number_unique = vec_unique(d_mesh->mesh_map(d_key)).size();
      Assert(number_unique == number_regions);

      // reshape the spectrum from 1-d to 2-d
      d_spectrum.resize(d_number_groups, vec_dbl(number_regions, 0.0));
      for (size_t g = 0; g < d_number_groups; ++g)
      {
        for (size_t r = 0; r < number_regions; ++r)
        {
          d_spectrum[g][r] = spectrum_1d[r + g * number_regions];
        }
      }
    }
    else if (d_condensation_option == CONDENSE_WITH_GS_SPECTRUM)
    {
      d_key = "MATERIAL";
      number_regions = d_material->number_materials();
      SpectrumGS S(d_input, d_material, d_mesh);
      d_spectrum = S.spectrum();
    }
    else if (d_condensation_option == CONDENSE_WITH_FS_SPECTRUM)
    {
      d_key = "MATERIAL";
      number_regions = d_material->number_materials();
      SpectrumFS S(d_input, d_material, d_mesh);
      d_spectrum = S.spectrum();
    }
    else if (d_condensation_option == CONDENSE_WITH_PINCELL_SPECTRUM)
    {
      d_key = "PINCELL_SPECTRUM";
      number_regions = d_material->number_materials();
      SpectrumPinCell S(d_input, d_material, d_mesh);
      d_spectrum = S.spectrum();
    }
    else
    {
      THROW("Invalid condensation option: " + AsString(d_condensation_option));
    }

    // Condense to coarse spectrum for shapes
    d_coarse_spectrum.resize(d_fine_per_coarse.size(),
                             vec_dbl(number_regions, 0.0));
    for (size_t r = 0; r < d_spectrum[0].size(); ++r)
    {
      size_t g_f = 0 ? d_adjoint : d_group_cutoff;
      for (size_t g_c = 0; g_c < d_fine_per_coarse.size(); ++g_c)
      {
        for (size_t i = 0; i < d_fine_per_coarse[g_c]; ++i, ++g_f)
        {
          d_coarse_spectrum[g_c][r] += d_spectrum[g_f][r];
        }
      }
    }
  }

  //--------------------------------------------------------------------------//
  // COARSE MATERIAL
  //--------------------------------------------------------------------------//

  Homogenize H(d_material);
  if (d_condensation_option == CONDENSE_WITH_STATE)
  {
    d_c_material = H.homogenize(state, d_mesh, "COARSEMESH",
                                d_fine_per_coarse, d_groups);
  }
  else
  {
    d_c_material = H.homogenize(d_spectrum, d_key, d_mesh, "COARSEMESH",
                                d_fine_per_coarse, d_groups);
//    std::cout << " SPECTRUM = " << std::endl;
//    for (int i = 0; i < d_spectrum.size(); ++i)
//    {
//      for (int j = 0; j < d_spectrum[i].size(); ++j)
//      {
//        std::cout << d_spectrum[i][j] << " ";
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    d_c_material->display();
  }

  // Prolongation
  build_prolong();

  }
  catch (detran_utilities::GenException &e)
  {
    std::cout << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << "Unknown exception!" << std::endl;
  }

  if (d_input->check("mgpc_print_operators"))
  {
    d_restrict->print_matlab("restrict.out");
    d_prolong->print_matlab("prolong.out");
  }

}

//----------------------------------------------------------------------------//
/*
 * Consider 3 groups with [2, 1] and 3 mesh with [2, 1].  The
 * restriction is then
 *
 *             R                  x  phi_f  =   phi_c
 *
 * | v0 v1  0 v0 v1  0  0  0  0 |   |f0g0|     |c0G0|
 * |  0  0 v3  0  0 v2  0  0  0 |   |f1g0|     |c1G0|
 * |  0  0  0  0  0  0 v0 v1  0 | x |f2g0|  =  |c0G1|
 * |  0  0  0  0  0  0  0  0 v2 |   |f0g1|     |c1G1|
 *                                  |f1g1|
 *                                  |f2g1|
 *                                  |f0g2|
 *                                  |f1g2|
 *                                  |f2g2|
 * where vk = v_i / v_j, where fine mesh i is in coarse mesh j.  Hence,
 * this produces a spatial average.  In energy, we simply add.  The two
 * combined preserve reaction rates as closely as possible.
 *
 * e.g   f = 1, g = 1  --> col = 3 * 1 + 1 = 4
 *                     --> row = 0 * 2 + 0 = 0
 */
void MGCoarseMeshPreconditioner::build_restrict()
{
  // Maximum number of nonzeros in a row.
  size_t nnz = (d_mesh->number_cells() / d_coarsemesh->number_cells() + 1) *
               detran_utilities::vec_max(d_fine_per_coarse);
  d_restrict = new callow::Matrix(d_size_coarse, d_size_fine, nnz);
  const vec_int &coarse_map = d_mesh->mesh_map("COARSEMESH");
  // Loop over active groups
  for (size_t g = 0; g < d_number_active_groups; ++g)
  {
    // Loop over fine mesh
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t c = coarse_map[i];
      size_t col = g * d_mesh->number_cells() + i;
      size_t row = d_f2c_group_map[g] * d_coarsemesh->number_cells() + c;
      double value = d_mesh->volume(i) / d_coarsemesh->volume(c);
      bool err = d_restrict->insert(row, col, value, callow::Matrix::INSERT);
      Assert(err);
    }
  }
  d_restrict->assemble();
}

//----------------------------------------------------------------------------//
/*
 * The prolongation maps the coarse fluxes back to their original
 * fine mesh/group location.  The simplest scheme is to set
 *    P = R'
 * and then set
 *   for all i, j do
 *     if (P(i, j) > 0) then P(i, j) = 1.0;
 * An alternative would be to use the spectral shape used to collapse
 * the groups.  For our 3 group example, suppose we have one representative
 * spectrum, [s0, s1, s2].  Mapping back from [2, 1] to [1, 1, 1] could
 * use [g0, g1, g2] <-- [s0/(s0+s1)*G0 s1/(s0+s1)*G0 s2/s2*G1]
 */
void MGCoarseMeshPreconditioner::build_prolong()
{
  // A fine mesh gets from only one coarse mesh
  d_prolong = new callow::Matrix(d_size_fine, d_size_coarse, 1);
  const vec_int &coarse_map = d_mesh->mesh_map("COARSEMESH");

  // Loop over active groups
  for (size_t g_f = 0; g_f < d_number_active_groups; ++g_f)
  {
    size_t g_c = d_f2c_group_map[g_f];

    // Loop over fine mesh
    for (size_t i_f = 0; i_f < d_mesh->number_cells(); ++i_f)
    {
      size_t i_c = coarse_map[i_f];
      size_t row = i_f + g_f * d_mesh->number_cells();
      size_t col = i_c + g_c * d_coarsemesh->number_cells();
      double value = shape(i_f, g_f, i_c, g_c);
      d_prolong->insert(row, col, value, callow::Matrix::INSERT);
    }
  }
  d_prolong->assemble();
}

//----------------------------------------------------------------------------//
double MGCoarseMeshPreconditioner::shape(const size_t i_f,
                                         const size_t g_f,
                                         const size_t i_c,
                                         const size_t g_c)
{
  // get the actual group to access fine group spectrum
  size_t g = d_adjoint ? g_f : g_f + d_group_cutoff;

  double value = 0.0;
  if (d_condensation_option == CONDENSE_WITH_STATE)
  {
    value = d_state->phi(g)[i_f] / d_coarse_spectrum[g_c][i_c];
  }
  else
  {
    const vec_int &spectrum_map = d_mesh->mesh_map(d_key);
    size_t i_s = spectrum_map[i_f];
    value = d_spectrum[g][i_s] / d_coarse_spectrum[g_c][i_s];
  }
  return value;
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGCoarseMeshPreconditioner.cc
//----------------------------------------------------------------------------//
