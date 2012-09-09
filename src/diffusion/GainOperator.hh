//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GainOperator.hh
 * \brief  GainOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef GAINOPERATOR_HH_
#define GAINOPERATOR_HH_

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC

// Diffusion
#include "BaseOperator.hh"
#include "utilities/Definitions.hh"

namespace detran_diffusion
{

/*!
 *  \class GainOperator
 *  \brief Gain operator for multigroup diffusion equation.
 *
 *  This operator is used in solution of multigroup diffusion
 *  problems.
 *
 *  A mesh-centered discretization is used.
 */
class GainOperator: public BaseOperator
{

public:

  typedef detran_utilities::SP<GainOperator>      SP_gainoperator;
  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_material::Material::SP_material  SP_material;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran_utilities::vec_int               vec_int;
  typedef detran_utilities::vec_dbl               vec_dbl;
  typedef detran_utilities::vec2_dbl              vec2_dbl;

  /*!
   *  \brief Constructor
   *  \param
   */
  GainOperator(SP_input    input,
               SP_material material,
               SP_mesh     mesh);

private:

  /// \name Private Data
  /// \{

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

#endif // GAINOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file GainOperator.hh
//---------------------------------------------------------------------------//
