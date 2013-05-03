//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   ScatterSource.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  ScatterSource class definition.
 */
//----------------------------------------------------------------------------//

#include "ScatterSource.hh"

namespace detran
{

//----------------------------------------------------------------------------//
ScatterSource::ScatterSource(SP_mesh     mesh,
                             SP_material material,
                             SP_state    state)
  :  d_mesh(mesh)
  ,  d_material(material)
  ,  d_state(state)
  ,  d_adjoint(false)
{
  Require(d_mesh);
  Require(d_material);
  Require(d_state);
  d_mat_map = d_mesh->mesh_map("MATERIAL");
  if (d_state->get_input()->check("adjoint"))
    d_adjoint = 0 != d_state->get_input()->get<int>("adjoint");
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of ScatterSource.cc
//----------------------------------------------------------------------------//
