//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Execute.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  Execute class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
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

// Config
#include "detran_config.h"

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "Eigensolver.hh"
#ifdef DETRAN_ENABLE_SLEPC
#include "EigenSLEPc.hh"
#endif
#include "ExternalSource.hh"
#include "FissionSource.hh"
#include "State.hh"
#include "StupidParser.hh"
#include "Traits.hh"
#include "PowerIteration.hh"
#include "GaussSeidel.hh"
#ifdef DETRAN_ENABLE_PETSC
#include "KrylovMG.hh"
#endif

// Utilities
#include "DBC.hh"
#include "InputDB.hh"

// System
#include <iostream>
#include <string>
#ifdef DETRAN_ENABLE_PETSC
#include "petsc.h"
#endif

namespace detran
{

//===========================================================================//
/*!
 * \class Execute
 * \brief Setup and execute the problem.
 */
//===========================================================================//

class Execute: public Object
{

public:

  // Basic Objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;


  // source typedefs
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;

  /*
   *  \brief Constructor
   *  \param parser     Input parser
   */
  Execute(StupidParser &parser);

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

  /// Unimplemented DBC method.
  bool is_valid() const
  {
    return true;
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
    boundary = new Boundary<D>(d_input, d_mesh, d_quadrature);

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
