//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGPreconditioner.cc
 *  @brief MGPreconditioner member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#include "MGPreconditioner.hh"
#include <cmath>

namespace detran
{

//----------------------------------------------------------------------------//
MGPreconditioner::MGPreconditioner(SP_input         input,
                                   SP_material      material,
                                   SP_mesh          mesh,
                                   SP_scattersource ssource,
                                   SP_fissionsource fsource,
                                   size_t           cutoff,
                                   bool             include_fission,
                                   bool             adjoint,
                                   std::string      name)
  : Base(name)
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_scattersource(ssource)
  , d_fissionsource(fsource)
  , d_group_cutoff(cutoff)
  , d_include_fission(include_fission)
  , d_adjoint(adjoint)
  , d_single_build(false)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  // The scatter and fission sources aren't necessarily required.

  // Number of groups
  d_number_groups = d_material->number_groups();

  // Number of active groups
  int upper = d_number_groups;
  if (d_adjoint) upper = -1;
  d_number_active_groups = std::abs(upper - (int)d_group_cutoff);

  // Check for a single build
  if (d_input->check("outer_pc_single_build"))
  {
    d_single_build = 0 != d_input->get<int>("outer_pc_single_build");
  }
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGPreconditioner.cc
//----------------------------------------------------------------------------//
