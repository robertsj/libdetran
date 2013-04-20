//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   WGPreconditioner.cc
 *  @brief  WGPreconditioner
 *  @author Jeremy Roberts
 *  @date   Nov 11, 2012
 */
//---------------------------------------------------------------------------//

#include "WGPreconditioner.hh"

namespace detran
{

//---------------------------------------------------------------------------//
WGPreconditioner::WGPreconditioner(SP_input input,
                                   SP_material material,
                                   SP_mesh mesh,
                                   std::string name)
  : Base(name)
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Number of groups
  d_number_groups = d_material->number_groups();

  // Size the solver and operator vectors
  d_solver.resize(d_number_groups);
  d_operator.resize(d_number_groups);

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file WGPreconditioner.cc
//---------------------------------------------------------------------------//
