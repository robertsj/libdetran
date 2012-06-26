//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryMOC.cc
 * \brief  BoundaryMOC member definitions.
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "BoundaryMOC.hh"
#include "VacuumMOC.hh"
#include "ReflectiveMOC.hh"

namespace detran
{

// Constructor.
template<class D>
BoundaryMOC<D>::BoundaryMOC(SP_input        input,
                            SP_mesh         mesh,
                            SP_quadrature   quadrature)
  : Base(input, mesh, quadrature)
  , d_boundary_flux(d_number_groups,
                    vec2_dbl(quadrature->number_angles()))
  , d_bc(2*D::dimension)
{
  // Allocated the flux container.
  initialize();

  // Create boundary conditions.
  std::vector<std::string> names(6);
  names[Mesh::LEFT]   = "bc_left";
  names[Mesh::RIGHT]  = "bc_right";
  names[Mesh::BOTTOM] = "bc_bottom";
  names[Mesh::TOP]    = "bc_top";
  names[Mesh::SOUTH]  = "bc_south";
  names[Mesh::NORTH]  = "bc_north";

  // Assign boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    // Vacuum is default.
    std::string type = "vacuum";
    if (input->check(names[side]))
    {
      type = input->get<std::string>(names[side]);
    }
    if (type == "vacuum")
    {
      d_bc[side] = new VacuumMOC<D>((*this), side, input, mesh, quadrature);
    }
    else if (type == "reflect")
    {
      d_bc[side] = new ReflectiveMOC<D>((*this), side, input, mesh, quadrature);
      d_is_reflective[side] = true;
      d_has_reflective = true;
    }
    else
    {
      type.append(" is not a supported bc type.");
      THROW(type);
      break;
    }
  }

}

template<class D>
void BoundaryMOC<D>::initialize()
{

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.cc
//---------------------------------------------------------------------------//
