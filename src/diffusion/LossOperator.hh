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

// Detran
#include "BaseOperator.hh"

// Utilities
#include "Definitions.hh"

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

  typedef detran::SP<LossOperator>      SP_lossoperator;
  typedef detran::InputDB::SP_input     SP_input;
  typedef detran::Material::SP_material SP_material;
  typedef detran::Mesh::SP_mesh         SP_mesh;
  typedef detran::vec_dbl               vec_dbl;
  typedef detran::vec2_dbl              vec2_dbl;

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

#endif // LOSSOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file LossOperator.hh
//---------------------------------------------------------------------------//
