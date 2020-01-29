//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ScatterSource.cc
 *  @brief ScatterSource class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
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
  d_adjoint = d_state->adjoint();
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of ScatterSource.cc
//----------------------------------------------------------------------------//
