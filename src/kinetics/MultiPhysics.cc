//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MultiPhysics.cc
 *  @brief  MultiPhysics
 *  @author Jeremy Roberts
 *  @date   Dec 19, 2012
 */
//---------------------------------------------------------------------------//

#include "MultiPhysics.hh"
#include "ioutils/StdOutUtils.hh"

namespace detran
{

//---------------------------------------------------------------------------//
MultiPhysics::MultiPhysics(size_t number_variables)
  : d_number_variables(number_variables)
  , d_physics_variables(number_variables)
{
  /* ... */
}

//---------------------------------------------------------------------------//
MultiPhysics::SP_multiphysics
MultiPhysics::Create(size_t number_variables)
{
  SP_multiphysics p(new MultiPhysics(number_variables));
  return p;
}

//---------------------------------------------------------------------------//
const MultiPhysics::vec_dbl&
MultiPhysics::variable(const size_t id) const
{
  // Preconditions
  Require(id < d_physics_variables.size());

  return d_physics_variables[id];
}

//---------------------------------------------------------------------------//
MultiPhysics::vec_dbl&
MultiPhysics::variable(const size_t id)
{
  // Preconditions
  Require(id < d_physics_variables.size());

  return d_physics_variables[id];
}

//---------------------------------------------------------------------------//
void MultiPhysics::add_variable(const size_t id, vec_dbl value)
{
  // Preconditions
  Require(id < d_physics_variables.size());

  d_physics_variables[id] = value;
}

//---------------------------------------------------------------------------//
void MultiPhysics::display() const
{
  /* ... */
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file MultiPhysics.cc
//---------------------------------------------------------------------------//
