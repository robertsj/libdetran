//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  WGDiffusionLossOperator.hh
 *  @brief WGDiffusionLossOperator class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_WGDIFFUSIONLOSSOPERATOR_HH_
#define detran_WGDIFFUSIONLOSSOPERATOR_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include "callow/matrix/Matrix.hh"

namespace detran
{

/**
 *  @class WGDiffusionLossOperator.hh
 *  @brief Loss operator for a one group diffusion equation.
 *
 *  This operator can be used for group sweeping in multigroup
 *  Krylov diffusion solvers.  Additionally, it can be used
 *  to precondition within-group transport solves.
 *
 *  The one group diffusion equation is
 *  \f[
 *      -\nabla D(\vec{r}) \nabla \phi(\vec{r})
 *        + \Sigma_r(\vec{r}) \phi(\vec{r}) = Q(\vec{r})
 *  \f]
 *  where the loss operator is defined
 *  \f[
 *    \mathbf{M}[\cdot] \equiv (-\nabla D(\vec{r}) \nabla
 *      + \Sigma_r(\vec{r}))[\cdot]\, .
 *  \f]
 *
 *  A mesh-centered discretization is used.
 */
class WGDiffusionLossOperator: public callow::Matrix
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef callow::Matrix                                Base;
  typedef callow::MatrixBase::SP_matrix                 SP_matrix;
  typedef detran_utilities::SP<WGDiffusionLossOperator> SP_operator;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_geometry::Mesh                         Mesh;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;
  typedef detran_utilities::size_t                      size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input            Pointer to input parameters
   *  @param material         Pointer to materials
   *  @param mesh             Pointer to mesh
   *  @param group            Group of operator
   */
  WGDiffusionLossOperator(SP_input    input,
                          SP_material material,
                          SP_mesh     mesh,
                          size_t      group);

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  /// Rebuild the matrix based on the present material definitions.
  void construct();

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
  /// Group of this operator
  size_t d_group;
  /// Problem dimension
  size_t d_dimension;
  /// Boundary albedo
  vec_dbl d_albedo;

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  // All matrix operator need this.  Here, we have the client call construct.
  void build();

};

} // end namespace detran

#endif // detran_WGDIFFUSIONLOSSOPERATOR_HH_

//---------------------------------------------------------------------------//
//              end of file WGDiffusionLossOperator.hh
//---------------------------------------------------------------------------//
