//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionSource.cc
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  FissionSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "FissionSource.hh"

// Utilities
#include "MathUtilities.hh"

// System
#include <iostream>

namespace detran
{

// Constructor
FissionSource::FissionSource(SP_state state,
                             SP_mesh mesh,
                             SP_material material)
  : d_state(state)
  , d_mesh(mesh)
  , d_material(material)
  , d_scale(1.0)
{
  Require(d_state);
  Require(d_mesh);
  Require(d_material);
  d_number_groups = d_material->number_groups();
  d_density.assign(d_mesh->number_cells(), 0.0);
  d_source.assign(d_mesh->number_cells(), 0.0);
}

void FissionSource::initialize()
{
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
  {
    d_density[cell] = d_material->nu_sigma_f(mat_map[cell], 0);
  }
  double norm_density = norm(d_density, "L2");
  Require(norm_density > 0.0);
  std::cout << " norm density = " << norm_density << std::endl;
  vec_scale(d_density, 1.0/norm_density);
}

void FissionSource::update()
{
  d_density.assign(d_density.size(), 0.0);
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for(int g = 0; g < d_number_groups; g++)
  {
    State::moments_type phi = d_state->phi(g);
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      d_density[cell] += phi[cell] *
                         d_material->nu_sigma_f(mat_map[cell], g);
    }
  }
}

void FissionSource::setup_outer(double scale)
{
  d_scale = scale;
}

const State::moments_type& FissionSource::source(int g)
{
  Require(g >= 0);
  Require(g < d_number_groups);
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
  {
    d_source[cell] = d_scale * d_density[cell] *
                     d_material->chi(mat_map[cell], g);
  }
  return d_source;
}

const State::moments_type& FissionSource::density()
{
  return d_density;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of FissionSource.cc
//---------------------------------------------------------------------------//
