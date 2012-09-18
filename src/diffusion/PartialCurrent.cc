//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PartialCurrent.cc
 * \author robertsj
 * \date   Sep 5, 2012
 * \brief  PartialCurrent class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "PartialCurrent.hh"
#include "DBC.hh"

namespace detran_diffusion
{

PartialCurrent::PartialCurrent(SP_input    input,
                               SP_material material,
                               SP_mesh     mesh)
 : d_input(input)
 , d_material(material)
 , d_mesh(mesh)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
}

void PartialCurrent::compute(Vec x, Vec b)
{

  // loop over groups
  for (int g = 0; g < d_material->number_groups(); g++)
  {





  } // end loop over groups

}


} // namespace detran_diffusion


