//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionGainOperator.hh
 * \brief  DiffusionGainOperator 
 * \author Jeremy Roberts
 * \date   Sep 10, 2012
 */
//---------------------------------------------------------------------------//

#ifndef DIFFUSIONGAINOPERATOR_HH_
#define DIFFUSIONGAINOPERATOR_HH_

#include "OperatorMatrix.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"

namespace detran
{

class DiffusionGainOperator: public OperatorMatrix
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<DiffusionGainOperator>   SP_gainoperator;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_geometry::Mesh                         Mesh;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param input            Pointer to input parameters
   *  \param material         Pointer to materials
   *  \param mesh             Pointer to mesh
   *  \param include_fission  Flag for including fission source implicitly
   *  \param keff             Fission scaling factor
   */
  DiffusionGainOperator(SP_input      input,
                        SP_material   material,
                        SP_mesh       mesh);

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

private:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

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

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  // All matrix operator need this.  Here, we have the client call construct.
  void build();

};

} // end namespace detran

#endif // DIFFUSIONGAINOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file DiffusionGainOperator.hh
//---------------------------------------------------------------------------//
