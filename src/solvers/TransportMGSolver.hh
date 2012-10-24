//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGTransportSolver.hh
 *  @author robertsj
 *  @date   Oct 24, 2012
 *  @brief  TransportMGSolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGTransportSolver_HH_
#define detran_MGTransportSolver_HH_

#include "InnerIteration.hh"
#include "angle/Quadrature.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/FissionSource.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  \class MGTransportSolver
 *  \brief Base class for multigroup transport solvers.
 *
 *  The multigroup transport equations are
 *  tbf.
 *
 */
//---------------------------------------------------------------------------//

template <class D>
class MGTransportSolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<SolverMG<D> >          SP_solver;
  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef detran_angle::Quadrature::SP_quadrature     SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_externalsource;
  typedef FissionSource::SP_fissionsource             SP_fissionsource;
  typedef typename InnerIteration<D>::SP_inner        SP_inner;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param external_source   User-defined external source.
   *  \param fission_source    Fission source.
   */
  SolverMG(SP_input           input,
           SP_state           state,
           SP_mesh            mesh,
           SP_material        material,
           SP_quadrature      quadrature,
           SP_boundary        boundary,
           SP_externalsource  q_e,
           SP_fissionsource   q_f);

  /// Virtual destructor
  virtual ~SolverMG(){};

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// Solve the multigroup equations.
  virtual void solve(const double keff = 1.0) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// User input.
  SP_input d_input;
  /// State vectors.
  SP_state d_state;
  /// Problem mesh
  SP_mesh d_mesh;
  /// Materials definitions.
  SP_material d_material;
  /// Angular mesh.
  SP_quadrature d_quadrature;
  /// Boundary fluxes
  SP_boundary d_boundary;
  /// External source
  SP_externalsource d_external_source;
  /// Fission source, if used
  SP_fissionsource d_fissionsource;
  /// Downscatter switch.
  bool d_downscatter;
  /// Number of energy groups.
  int d_number_groups;
  /// Maximum outer iterations (only relevant for upscatter)
  int d_max_iters;
  /// Outer tolerance
  double d_tolerance;
  /// Print out flag
  int d_print_out;
  /// Interval for print out
  int d_print_interval;
  /// Inner solver
  SP_inner d_inner_solver;
  /// Multiplying problem flag
  bool d_multiply;

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "SolverMG.i.hh"


#endif /* detran_TRANSPORTMGSOLVER_HH_ */
