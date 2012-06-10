//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ScatterSource.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  ScatterSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "ScatterSource.hh"

namespace detran
{

ScatterSource::ScatterSource(SP_mesh mesh,
                             SP_material material,
                             SP_state state)
  :  d_mesh(mesh)
  ,  d_material(material)
  ,  d_state(state)
{
  Require(d_mesh);
  Require(d_material);
  Require(d_state);

  // \todo Add a check function to mesh like input has.
  d_mat_map = d_mesh->mesh_map("MATERIAL");

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of ScatterSource.cc
//---------------------------------------------------------------------------//
