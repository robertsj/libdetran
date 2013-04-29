//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryFactory.hh
 *  @author robertsj
 *  @date   Jan 30, 2013
 *  @brief  BoundaryFactory class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYFACTORY_HH_
#define detran_BOUNDARYFACTORY_HH_

#include "utilities/InputDB.hh"
#include "geometry/Mesh.hh"
#include "angle/Quadrature.hh"

namespace detran
{

/// Class for constructing boundaries and boundary conditions
template <typename D, template <typename> class B>
class BoundaryFactory
{
    
public:
    
  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef typename B<D>::SP_boundary                    SP_boundary;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature       SP_quadrature;
    
  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  static SP_boundary build(SP_input input, SP_mesh mesh, SP_quadrature quad)
  {
    THROW("NOT IMPLEMENTED");
    return SP_boundary();
  }
};

} // end namespace detran

#endif /* detran_BOUNDARYFACTORY_HH_ */
