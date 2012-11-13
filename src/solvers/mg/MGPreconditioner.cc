//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGPreconditioner.cc
 *  @brief  MGPreconditioner
 *  @author Jeremy Roberts
 *  @date   Nov 12, 2012
 */
//---------------------------------------------------------------------------//


#include "MGPreconditioner.hh"

namespace detran
{

//---------------------------------------------------------------------------//
MGPreconditioner::MGPreconditioner(SP_input input,
                                   SP_material material,
                                   SP_mesh mesh,
                                   size_t cutoff,
                                   std::string name)
  : Base(name)
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_group_cutoff(cutoff)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Number of groups
  d_number_groups = d_material->number_groups();

  // Number of active groups
  d_number_active_groups = d_number_groups - d_group_cutoff;

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file MGPreconditioner.cc
//---------------------------------------------------------------------------//
