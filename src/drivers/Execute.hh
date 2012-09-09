//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Execute.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  Execute class definition.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! \mainpage detran: A DEterministic TRANsport package
 *
 * \section sec_introduction Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */
//---------------------------------------------------------------------------//

#ifndef EXECUTE_HH_
#define EXECUTE_HH_

#include "detran_config.h"
#include "StupidParser.hh"
#include "discretization/DimensionTraits.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "solvers/GaussSeidel.hh"
#include "solvers/Eigensolver.hh"
#ifdef DETRAN_ENABLE_SLEPC
#include "solvers/EigenSLEPc.hh"
#endif
#ifdef DETRAN_ENABLE_PETSC
#include "solvers/KrylovMG.hh"
#include "petsc.h"
#endif
#include "solvers/PowerIteration.hh"
#include "transport/FissionSource.hh"
#include "transport/State.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include <iostream>
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Execute
 * \brief Setup and execute the problem.
 */
//---------------------------------------------------------------------------//

class Execute
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input             SP_input;
  typedef State::SP_state                                 SP_state;
  typedef detran_geometry::Mesh::SP_mesh                  SP_mesh;
  typedef detran_material::Material::SP_material          SP_material;
  typedef detran_angle::Quadrature::SP_quadrature         SP_quadrature;
  typedef detran_external_source::
          ExternalSource::SP_externalsource               SP_externalsource;
  typedef detran::FissionSource::SP_fissionsource         SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*
   *  \brief Constructor
   *  \param parser     Input parser
   */
  Execute(StupidParser &parser);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Solve the problem.
  template <class D>
  void solve();

  /// Write output to file.
  void output();

  // Return dimension of problem.
  int dimension()
  {
    return d_dimension;
  }

private:

  // Input
  SP_input d_input;
  // Material
  SP_material d_material;
  // Mesh
  SP_mesh d_mesh;
  // State
  SP_state d_state;
  //
  SP_quadrature d_quadrature;
  //
  SP_externalsource d_externalsource;
  //
  SP_fissionsource d_fissionsource;
  //
  int d_dimension;
  //
  std::string d_problem_type;
  //
  bool d_moc;

  /// Do nontemplated portion of problem setup.
  void setup();

};

template<class D>
void Execute::solve()
{

  //--------------------------------------------------------------------------//
  // Boundary
  //--------------------------------------------------------------------------//

  typename BoundaryBase<D>::SP_boundary boundary;
  if (d_moc)
    boundary = new BoundaryMOC<D>(d_input, d_mesh, d_quadrature);
  else
    boundary = new BoundarySN<D>(d_input, d_mesh, d_quadrature);

  //--------------------------------------------------------------------------//
  // Create solver and solve.
  //--------------------------------------------------------------------------//

  if (d_problem_type == "eigenvalue")
  {
    std::string eigen_solver = "PI";
    if (d_input->check("eigen_solver"))
    {
      eigen_solver = d_input->get<std::string>("eigen_solver");
    }

    if (eigen_solver == "PI")
    {
      PowerIteration<D> solver(d_input, d_state, d_mesh, d_material,
                               d_quadrature, boundary, d_fissionsource);
      solver.solve();
    }
    else if (eigen_solver == "SLEPc")
    {
#ifdef DETRAN_ENABLE_SLEPC
      EigenSLEPc<D> solver(d_input, d_state, d_mesh, d_material,
                           d_quadrature, boundary, d_fissionsource);
      solver.solve();
#else
      THROW("EigenSLEPc is unavailable since SLEPc is not enabled.");
#endif
    }
    else
    {
      THROW("Unsupported eigen_solver type selected: "+eigen_solver);
    }

  }
  else if (d_problem_type == "fixed" || d_problem_type == "fixed_multiply")
  {
    std::string outer_solver = "GS";
    if (d_input->check("outer_solver"))
    {
      outer_solver = d_input->get<std::string>("outer_solver");
    }
    if (outer_solver == "GS")
    {
      GaussSeidel<D> solver(d_input, d_state, d_mesh, d_material, d_quadrature,
                            boundary, d_externalsource, d_fissionsource);
      solver.solve();
    }
    else if (outer_solver == "KrylovMG")
    {
#ifdef DETRAN_ENABLE_PETSC
      KrylovMG<D> solver(d_input, d_state, d_mesh, d_material, d_quadrature,
                         boundary, d_externalsource, d_fissionsource);
      solver.solve();
#else
      THROW("KrylovMG is unavailable since PETSc is not enabled.");
#endif
    }
    else
    {
      THROW("Unsupported outer_solver type selected: "+outer_solver);
    }
  }
  else
  {
    std::cout << "Unsupported problem type given.  Options are:" << std::endl
              << " -- fixed" << std::endl
              << " -- eigenvalue" << std::endl
              << std::endl;
  }

  // INSERT POST PROCESSING HERE
  State::moments_type phi = d_state->phi(0);
  std::cout << phi[0] << " " << phi[1] << std::endl;

}

} // end namespace detran

#endif /* EXECUTE_HH_ */

//---------------------------------------------------------------------------//
//              end of Execute.hh
//---------------------------------------------------------------------------//
