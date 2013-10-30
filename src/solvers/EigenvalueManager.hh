//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenvalueManager.hh
 *  @brief EigenvalueManager class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENVALUEMANAGER_HH_
#define detran_EIGENVALUEMANAGER_HH_

#include "TransportManager.hh"
#include "solvers/eigen/Eigensolver.hh"

namespace detran
{

/**
 *  @class EigenvalueManager
 *  @brief Manage solution of a multigroup eigenvalue problem
 */
template <class D>
class EigenvalueManager: public TransportManager
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef typename Eigensolver<D>::Fixed_T              Fixed_T;
  typedef typename Eigensolver<D>::SP_solver            SP_solver;
  typedef typename Fixed_T::SP_manager                  SP_mg_solver;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef State::SP_state                               SP_state;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_angle::Quadrature::SP_quadrature       SP_quadrature;
  typedef BoundaryBase<D>                               Boundary_T;
  typedef typename Boundary_T::SP_boundary              SP_boundary;
  typedef FissionSource::SP_fissionsource               SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param argc       command line count
   *  @param argv       command line values
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   */
  EigenvalueManager(int argc,
                    char *argv[],
                    SP_input    input,
                    SP_material material,
                    SP_mesh     mesh);

  /// Constructor (without command line)
  EigenvalueManager(SP_input    input,
                    SP_material material,
                    SP_mesh     mesh);

  /// Virtual destructor
  virtual ~EigenvalueManager(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Solve the system
  bool solve();

  /// Getters
  SP_input input() const {return d_mg_solver->input();}
  SP_material material() const {return d_mg_solver->material();}
  SP_mesh mesh() const {return d_mg_solver->mesh();}
  SP_state state() const {return d_mg_solver->state();}
  SP_boundary boundary() const {return d_mg_solver->boundary();}
  SP_quadrature quadrature() const {return d_mg_solver->quadrature();}
  SP_fissionsource fissionsource() const {return d_mg_solver->fissionsource();}
  int number_sweeps() const { return d_mg_solver->number_sweeps(); }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Eigensolver
  SP_solver d_solver;
  /// Multigroup solver
  SP_mg_solver d_mg_solver;
  /// Adjoint mode flag
  bool d_adjoint;
  /// Discretization general type (SN, MOC, or diffusion)
  int d_discretization;
  /// Problem setup status flag
  bool d_is_setup;

};

} // end namespace detran

#endif /* detran_EIGENVALUEMANAGER_HH_ */

//----------------------------------------------------------------------------//
//              end of EigenvalueManager.hh
//----------------------------------------------------------------------------//
