//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PyExecute.hh
 * \author robertsj
 * \date   Jun 13, 2012
 * \brief  PyExecute class definition.
 */
//---------------------------------------------------------------------------//

#ifndef PYEXECUTE_HH_
#define PYEXECUTE_HH_

// Config
#include "detran_config.h"

// Detran
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
 * \class PyExecute
 * \brief Setup and execute the problem from the Python front end.
 */
//---------------------------------------------------------------------------//

template <class D>
class PyExecute
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utility::InputDB::SP_input       SP_input;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran_material::Material::SP_material  SP_material;

  typedef State::SP_state                         SP_state;
  typedef detran_angle::Quadrature::SP_quadrature SP_quadrature;
  typedef detran_external_source::
          ExternalSource::SP_externalsource       SP_externalsource;
  typedef FissionSource::SP_fissionsource         SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*
   *  \brief Constructor
   *
   *  This initializes PETSc, if applicable.
   *
   *  \param argc   Command line argument count
   *  \param argv   Command line arguments
   */
  PyExecute(int argc, char *argv[]);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /*
   *  \brief Initialize the manager for a new transport problem.
   *
   *  If a new problem is to be solved in the same script, a user
   *  need only call initialize again, thus resetting all the
   *  internal pointers.
   *
   *  \param input    User input database
   *  \param material Cross section library
   *  \param mesh     Problem geometry
   */
  void initialize(SP_input input, SP_material material, SP_mesh mesh);

  /// Solve the problem.
  void solve();

  /// Finalize.  Closes PETSc, if running, etc.
  void finalize();

  /// \name Setters
  /// \{

  /// Set an external source.
  void set_external_source(SP_externalsource q)
  {
    Require(q);
    d_externalsource = q;
  }

  /// \}

  /// \name Getters
  /// \{

  SP_state get_state()
  {
    return d_state;
  }

  SP_quadrature get_quadrature()
  {
    return d_quadrature;
  }

  /// \}

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  SP_input d_input;
  SP_material d_material;
  SP_mesh d_mesh;
  SP_state d_state;
  SP_quadrature d_quadrature;
  SP_externalsource d_externalsource;
  SP_fissionsource d_fissionsource;
  std::string d_problem_type;
  bool d_initialized;
  bool d_moc;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Generic problem setup.
  void setup();

};

template <class D>
PyExecute<D>::PyExecute(int argc, char *argv[])
  : d_initialized(false)
{
#ifdef DETRAN_ENABLE_PETSC
  // Start PETSc.
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif
#ifdef DETRAN_ENABLE_SLEPC
  // Start SLEPc.
  SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif
}

template <class D>
void PyExecute<D>::initialize(SP_input input,
                              SP_material material,
                              SP_mesh mesh)
{
  // Check preconditions.
  Require(input);
  Require(material);
  Require(mesh);
  d_input = input;
  d_material = material;
  d_mesh = mesh;

  // Setup the problem.
  setup();
  d_initialized = true;

  // Check postconditions.
  Ensure(d_quadrature);
}

template <class D>
void PyExecute<D>::finalize()
{
#ifdef DETRAN_ENABLE_PETSC
  PetscFinalize();
#endif
#ifdef DETRAN_ENABLE_SLEPC
  SlepcFinalize(); // is this safe?
#endif
}

//---------------------------------------------------------------------------//
// Private implementation.
//---------------------------------------------------------------------------//

template <class D>
void PyExecute<D>::setup()
{
  using std::string;

  // Ensure material and input groups match.
  d_input->put<int>("number_groups", d_material->number_groups());

  // Problem type
  if (d_input->check("problem_type"))
    d_problem_type = d_input->get<std::string>("problem_type");

  // Equation type and MOC flag.
  string eq = "dd";
  d_moc = false;
  if (d_input->check("equation"))
    eq = d_input->get<std::string>("equation");
  if (eq == "scmoc" || eq == "ddmoc")
    d_moc = true;

  //-------------------------------------------------------------------------//
  // Quadrature
  //-------------------------------------------------------------------------//

  QuadratureFactory quad_factory;
  quad_factory.build(d_quadrature, d_input, D::dimension);
  Assert(d_quadrature);

  //-------------------------------------------------------------------------//
  // MOC Mesh
  //-------------------------------------------------------------------------//
  if (d_moc)
  {
    // Track the mesh
    Tracker tracker(d_mesh, d_quadrature);

    // Normalize segments to conserve volume.
    tracker.normalize();

    // Replace the mesh with the tracked one.  This suggests refactoring
    // to have a (possibly null) trackdb in Mesh.
    d_mesh = tracker.meshmoc();
  }

  //-------------------------------------------------------------------------//
  // State
  //-------------------------------------------------------------------------//

  d_state = new detran::State(d_input, d_mesh, d_quadrature);

  //-------------------------------------------------------------------------//
  // Fission source
  //-------------------------------------------------------------------------//

  if (d_problem_type == "eigenvalue" || d_problem_type == "fixed_multiply")
  {
    d_fissionsource = new FissionSource(d_state, d_mesh, d_material);
  }

}

template<class D>
void PyExecute<D>::solve()
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
    Require(d_fissionsource);

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

}

// Explicit instantiations
template class PyExecute<_1D>;
template class PyExecute<_2D>;
template class PyExecute<_3D>;

} // end namespace detran

#endif /* PYEXECUTE_HH_ */

//---------------------------------------------------------------------------//
//              end of PyExecute.hh
//---------------------------------------------------------------------------//
