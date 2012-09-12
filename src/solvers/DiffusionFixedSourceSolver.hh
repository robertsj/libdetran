//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionFixedSourceSolver.hh
 * \brief  DiffusionFixedSourceSolver class definition
 * \author Jeremy Roberts
 * \date   Sep 11, 2012
 */
//---------------------------------------------------------------------------//

#ifndef DIFFUSIONFIXEDSOURCESOLVER_HH_
#define DIFFUSIONFIXEDSOURCESOLVER_HH_

#include "DiffusionLossOperator.hh"
#include "DiffusionGainOperator.hh"
#include "Vector.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "external_source/IsotropicSource.hh"
#include "utilities/Definitions.hh"
#include "petsc.h"

namespace detran
{

/*!
 *  \class DiffusionFixedSourceSolver
 *  \brief Solve a fixed source problem using diffusion
 *
 *  Fixed source problems can be solved in the following modalities:
 *    - fixed, no multiplication
 *    - fixed, multiplying via implicit fission (i.e. in the loss matrix)
 *    - fixed, multiplying via fission iteration
 */
template <class D>
class DiffusionFixedSourceSolver
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum fixed_types
  {
    FIXED, MULTIPLY, ITERATE, END_FIXED_TYPES
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef typename BoundaryDiffusion<D>::SP_boundary  SP_boundary;
  typedef detran_external_source::
          IsotropicSource::SP_externalsource          SP_source;
  typedef State::moments_type                         moments_type;
  typedef DiffusionLossOperator::SP_lossoperator      SP_lossoperator;
  typedef DiffusionGainOperator::SP_gainoperator      SP_gainoperator;
  typedef Vector::SP_vector                           SP_vector;
  typedef detran_utilities::size_t                    size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  DiffusionFixedSourceSolver(SP_input input,
                             SP_material material,
                             SP_mesh mesh,
                             SP_state state);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  void solve();

  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  SP_input d_input;
  SP_material d_material;
  SP_mesh d_mesh;
  SP_state d_state;
  SP_boundary d_boundary;
  SP_lossoperator d_M;
  SP_gainoperator d_F;
  SP_vector d_x;
  SP_vector d_b;
  size_t d_problem_size;
  size_t d_fixed_type;
  size_t d_maximum_iterations;
  size_t d_maximum_fission_iterations;
  double d_tolerance;
  double d_fission_tolerance;
  double d_fission_scaling;
  KSP d_solver;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  void solve_fixed();
  void solve_iterate();

};

} // end namespace detran

#endif // DIFFUSIONFIXEDSOURCESOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file DiffusionFixedSourceSolver.hh
//---------------------------------------------------------------------------//
