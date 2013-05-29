//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCoarseMeshPreconditioner.cc
 *  @brief MGCoarseMeshPreconditioner member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGCoarseMeshPreconditioner.hh"
#include "transport/Homogenize.hh"
#include "utilities/MathUtilities.hh"

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
{
  /* ... */
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

  //--------------------------------------------------------------------------//
  // SPECTRUM FOR WEIGHTING
  //--------------------------------------------------------------------------//

  size_t option = CONDENSE_WITH_UNITY_SPECTRUM;
  if (d_input->check("mgpc_condensation_option"))
  {
    option = d_input->get<int>("mgpc_condensation_option");
    Insist(option < END_CONDENSATION_OPTIONS,
           "Invalid MGPC condensation option selected.");
  }
  vec2_dbl spectrum;
  std::string key;
  if (option == CONDENSE_WITH_STATE)
  {
    Insist(state, "To build a MGPC using a state requires state != NULL");
  }
  else
  {
    // Only unity weighting for now.  Note, the spectrum has to have all
    // groups, even though the material might not be homogenized over all
    // groups.
    spectrum.resize(d_number_groups,
                    vec_dbl(d_material->number_materials(), 1.0));
    key = "MATERIAL";
  }

  //--------------------------------------------------------------------------//
  // COARSE MESH
  //--------------------------------------------------------------------------//

  size_t level = 2;
  if (d_input->check("mgpc_coarseness_level"))
  {
    level = d_input->get<int>("mgpc_coarseness_level");
    Assert(level > 0);
  }
  CoarseMesh C(d_mesh, level);
  d_coarsemesh = C.get_coarse_mesh();

  //--------------------------------------------------------------------------//
  // COARSE MATERIAL
  //--------------------------------------------------------------------------//

  Homogenize H(d_material);

  // fine to coarse map --- defined *only* for active groups
  vec_int coarsegroups(d_number_active_groups, 1);

  int upper = d_adjoint ? -1 : d_number_groups;
  groups_t groups = detran_utilities::range<size_t>(d_group_cutoff, upper);

  if (d_input->check("mgpc_coarsegroups"))
  {
    coarsegroups = d_input->get<vec_int>("mgpc_coarsegroups");
    Require(detran_utilities::vec_sum(coarsegroups) == d_number_groups);
  }
  if (option == CONDENSE_WITH_STATE)
  {
    d_coarsematerial =
      H.homogenize(state, d_mesh, "COARSEMESH", coarsegroups, groups);
  }
  else
  {
    d_coarsematerial =
      H.homogenize(spectrum, key, d_mesh, "COARSEMESH", coarsegroups, groups);
  }

  //--------------------------------------------------------------------------//
  // RESTRICTION AND PROLONGATION
  //--------------------------------------------------------------------------//

  /*
   *  For restriction and prolongation, we need to account for the possible
   *  trucated groups.  Same for homogenization.
   */

  size_t size_fine   = d_number_active_groups * d_mesh->number_cells();
  size_t size_coarse = d_number_active_groups * d_coarsemesh->number_cells();

  // **** NEED TO INCLUDE ENERGY!!!

  // Restriction
  const vec_int &reg_map = d_mesh->mesh_map("COARSEMESH");
  d_restrict = new callow::Matrix(size_coarse, size_fine, level+1);
  for (size_t fine = 0; fine < size_fine; ++fine)
  {
    size_t coarse = reg_map[fine];
    double value = 1.0 / d_coarsemesh->volume(coarse);
    d_restrict->insert(fine, coarse, value, callow::Matrix::INSERT);
  }
  d_restrict->assemble();

  // Prolongation
  d_prolong = new callow::Matrix(size_fine, size_coarse, 1);
  for (size_t fine = 0; fine < size_fine; ++fine)
  {
    size_t coarse = reg_map[fine];
    d_restrict->insert(coarse, fine, 1.0, callow::Matrix::INSERT);
  }
  d_prolong->assemble();
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGCoarseMeshPreconditioner.cc
//----------------------------------------------------------------------------//
