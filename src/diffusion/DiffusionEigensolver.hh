//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionEigensolver.hh
 * \author robertsj
 * \date   Jul 26, 2012
 * \brief  DiffusionEigensolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef DIFFUSIONEIGENSOLVER_HH_
#define DIFFUSIONEIGENSOLVER_HH_

// Diffusion
#include "LossOperator.hh"
#include "GainOperator.hh"

// Transport
#include "State.hh"

// Utilities
#include "Definitions.hh"

// System
#include "slepc.h"

namespace detran_diffusion
{

/*!
 *  \class DiffusionEigensolver
 *  \brief Eigensolver for multigroup diffusion equations.
 */
class DiffusionEigensolver
{

public:

  typedef detran::InputDB::SP_input     SP_input;
  typedef detran::Material::SP_material SP_material;
  typedef detran::Mesh::SP_mesh         SP_mesh;
  typedef detran::State::SP_state       SP_state;
  typedef LossOperator::SP_lossoperator SP_lossoperator;
  typedef GainOperator::SP_gainoperator SP_gainoperator;

  /*!
   *  \brief Constructor
   *  \param
   */
  DiffusionEigensolver(SP_input    input,
                       SP_material material,
                       SP_mesh     mesh,
                       SP_state    state);

  /// Solve the system
  void solve();

  /// Get the loss operator
  detran::SP<detran_diffusion::LossOperator> M()
  {
    return d_M;
  }

  /// Get the gain operator
  detran::SP<detran_diffusion::GainOperator> F()
  {
    return d_F;
  }

private:

  /// \name Private Data
  /// \{

  /// Loss operator
  SP_lossoperator d_M;

  /// Gain operator
  SP_gainoperator d_F;

  /// SLEPc eigensolver
  EPS d_solver;

  /// Eigenvector
  Vec d_x;

  /// Eigenvalue
  double d_keff;

  /// System size
  int d_size;

  /// Number of groups
  int d_number_groups;

  /// Group size
  int d_group_size;

  /// Input database
  SP_input d_input;

  /// State
  SP_state d_state;

  /// Maximum eigensolver iterations
  int d_max_iters;

  /// Eigensolver tolerance, \f$ ||\mathbf{A} x-\lambda x|| < tol \f$.
  double d_tolerance;

  /// \}

  /// \name Implementation
  /// \{

  /// Fill the state with the eigenvector and the eigenvalue.
  void fill_state();

  /// \}

};

} // end namespace detran_diffusion

#endif /* DIFFUSIONEIGENSOLVER_HH_ */
