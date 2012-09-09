//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LossOperator.hh
 * \brief  LossOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef LOSSOPERATOR_HH_
#define LOSSOPERATOR_HH_

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC

// Diffusion
#include "BaseOperator.hh"

// Utilities
#include "utilities/Definitions.hh"

namespace detran_diffusion
{

/*!
 *  \class LossOperator
 *  \brief Loss operator for multigroup diffusion equation.
 *
 *  This operator is used in solution of multigroup diffusion
 *  problems.
 *
 *  A mesh-centered discretization is used.
 */
class LossOperator: public BaseOperator
{

public:

  typedef detran_utilities::SP<LossOperator>      SP_lossoperator;
  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_material::Material::SP_material  SP_material;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran_geometry::Mesh                   Mesh;
  typedef detran_utilities::vec_int               vec_int;
  typedef detran_utilities::vec_dbl               vec_dbl;
  typedef detran_utilities::vec2_dbl              vec2_dbl;

  /*!
   *  \brief Constructor
   *  \param
   */
  LossOperator(SP_input    input,
               SP_material material,
               SP_mesh     mesh);

private:

  /// \name Private Data
  /// \{

  /// Energy-dependent albedos
  vec2_dbl d_albedo;

  /// Energy groups
  int d_number_groups;

  /// One group spatial size
  int d_group_size;

  /// \}

  /// \name Implementation
  /// \{

  void construct();

  /// \}

};

} // end namespace detran_diffusion

#endif // DETRAN_ENABLE_PETSC

#endif // LOSSOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file LossOperator.hh
//---------------------------------------------------------------------------//
