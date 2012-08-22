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
#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC

// Detran
#include "Material.hh"
#include "Mesh.hh"

// Utilities
#include "InputDB.hh"
#include "Timer.hh"
#include "SP.hh"

// System
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

  typedef detran::SP<BaseOperator>      SP_operator;
  typedef detran::InputDB::SP_input     SP_input;
  typedef detran::Material::SP_material SP_material;
  typedef detran::Mesh::SP_mesh         SP_mesh;

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
  detran::Timer d_timer;

  /// \}

private:

};

} // end namespace detran

#endif // DETRAN_ENABLE_PETSC

#endif // BASEOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file BaseOperator.hh
//---------------------------------------------------------------------------//
