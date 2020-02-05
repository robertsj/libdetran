//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DiffusionGainOperator.hh
 *  @brief DiffusionGainOperator
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_DIFFUSIONGAINOPERATOR_HH_
#define detran_DIFFUSIONGAINOPERATOR_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include "callow/matrix/Matrix.hh"

namespace detran
{

//----------------------------------------------------------------------------//
class DiffusionGainOperator: public callow::Matrix
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef callow::Matrix                                Base;
  typedef detran_utilities::SP<DiffusionGainOperator>   SP_gainoperator;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_geometry::Mesh                         Mesh;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   */
  DiffusionGainOperator(SP_input      input,
                        SP_material   material,
                        SP_mesh       mesh,
                        bool          adjoint = false);

  /// SP constructor
  static SP_gainoperator Create(SP_input    input,
                                SP_material material,
                                SP_mesh     mesh,
                                bool        adjoint = false)
  {
    SP_gainoperator p(new DiffusionGainOperator(input, material, mesh, adjoint));
    return p;
  }

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  void construct(SP_material mat);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Input
  SP_input d_input;
  /// Material
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// Problem dimension
  size_t d_dimension;
  /// Energy groups
  size_t d_number_groups;
  /// One group spatial size
  size_t d_group_size;
  /// Adjoint flag
  bool d_adjoint;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  // All matrix operator need this.  Here, we have the client call construct.
  void build();

};

} // end namespace detran

#endif // detran_DIFFUSIONGAINOPERATOR_HH_

//----------------------------------------------------------------------------//
//              end of file DiffusionGainOperator.hh
//----------------------------------------------------------------------------//
