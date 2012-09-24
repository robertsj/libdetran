//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BaseOperator.hh
 * \brief  BaseOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BASEOPERATOR_HH_
#define BASEOPERATOR_HH_

// Configuration
#include "config/detran_config.hh"

#ifdef DETRAN_ENABLE_PETSC
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "utilities/Timer.hh"
#include "utilities/SP.hh"
#include "petsc.h"

namespace detran_diffusion
{

/*!
 *  \class BaseOperator
 *  \brief Base operator for a diffusion equation.
 *
 */
class BaseOperator
{

public:

  typedef detran_utilities::SP<BaseOperator>      SP_operator;
  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_material::Material::SP_material  SP_material;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;

  /*!
   *  \brief Constructor
   *
   */
  BaseOperator(SP_input    input,
               SP_material material,
               SP_mesh     mesh);

  /// Virtual destructor
  virtual ~BaseOperator();

  /// Return the PETSc operator
  Mat get_operator()
  {
    return d_operator;
  }

  /*!
   *  \brief Matrix-vector multiplication
   *  \param v_in   Input vector
   *  \param v_out  Output vector
   */
  void multiply(Vec v_in, Vec v_out);

  /*!
   *  \brief Print to standard out
   */
  void display();

protected:

  /// \name Protected Data
  /// \{

  /// PETSc matrix
  Mat d_operator;

  /// Size of the operator
  int d_size;

  /// Input
  SP_input d_input;

  /// Materials
  SP_material d_material;

  /// Mesh
  SP_mesh d_mesh;

  /// Dimension
  int d_dimension;

  /// Timer
  detran_utilities::Timer d_timer;

  /// \}

private:

};

} // end namespace detran

#endif // DETRAN_ENABLE_PETSC

#endif // BASEOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file BaseOperator.hh
//---------------------------------------------------------------------------//
