//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Solver.hh
 * \brief  Solver class definition
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef SOLVER_HH_
#define SOLVER_HH_

#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"
#include "geometry/Mesh.hh"

namespace detran
{

/*!
 *  \class Solver
 *  \brief Base class for transport solvers
 *
 *
 */
template <class D>
class Solver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran_material::Material::SP_material  SP_material;
  typedef State::SP_state                         SP_state;

  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;

  // source typedefs
  typedef ExternalSource::SP_source SP_externalsource;
  typedef FissionSource::SP_source SP_fissionsource;
  // sweep
  typedef typename Sweeper<D>::SP_sweeper SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource SP_sweepsource;
  // acceleration
  //  typedef typename
  //      Acceleration<D>::SP_acceleration          SP_acceleration;
  //
  typedef State::moments_type moments_type;


  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  Solver();

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SOLVERS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

protected:

  /// User input.
  SP_input d_input;
  /// State vectors.
  SP_state d_state;
  /// Problem mesh (either Cartesian mesh or MOC tracking)
  SP_mesh d_mesh;
  /// Materials definitions.
  SP_material d_material;

  /// Boundary fluxes.
  SP_boundary d_boundary;

  /// User-defined external source
  SP_externalsource d_external_source;
  /// Fission source, if used
  SP_fissionsource d_fission_source;

  /// Maximum iterations
  size_t d_max_iterations;

  /// Convergence tolerance
  double d_tolerance;

  /// Print out flag
  int d_print_out;
  /// Interval for print out
  int d_print_interval;

};

} // end namespace detran

#endif // SOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file Solver.hh
//---------------------------------------------------------------------------//
