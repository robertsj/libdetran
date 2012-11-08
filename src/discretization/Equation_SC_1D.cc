//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_1D.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Equation_SC_1D class definition.
 */
//---------------------------------------------------------------------------//

#include "Equation_SC_1D.hh"

namespace detran
{

//---------------------------------------------------------------------------//
Equation_SC_1D::Equation_SC_1D(SP_mesh mesh,
                               SP_material material,
                               SP_quadrature quadrature,
                               bool update_psi)
  :  Equation<_1D>(mesh, material, quadrature, update_psi)
  ,  d_mu(-1.0)
{
  /* ... */
}

//---------------------------------------------------------------------------//
void Equation_SC_1D::setup_group(const size_t g)
{
  Require(g < d_material->number_groups());
  d_g = g;
}

//---------------------------------------------------------------------------//
void Equation_SC_1D::setup_octant(const size_t octant)
{
  Require(octant < 2);
  d_octant = octant;
}

} // end namespace detran


