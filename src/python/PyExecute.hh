//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PyExecute.hh
 * \author robertsj
 * \date   Jun 13, 2012
 * \brief  PyExecute class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef PYEXECUTE_HH_
#define PYEXECUTE_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "ExternalSource.hh"
#include "FissionSource.hh"
#include "State.hh"
#include "StupidParser.hh"
#include "Traits.hh"
#include "PowerIteration.hh"
#include "GaussSeidel.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"

// System
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
class PyExecute: public Object
{

public:

  /// \name Useful Typedefs
  /// \{

  typedef InputDB::SP_input                     SP_input;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;

  typedef State::SP_state                       SP_state;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;

  /*
   *  \brief Constructor
   *  \param input    User input database
   *  \param material Cross section library
   *  \param mesh     Problem geometry
   */
  PyExecute(SP_input input, SP_material material, SP_mesh mesh);

  /// Solve the problem.
  void solve();

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

  /// Unimplemented DBC method.
  bool is_valid() const
  {
    return true;
  }

private:

  /// \name Private Data
  /// \{

  SP_input d_input;
  SP_material d_material;
  SP_mesh d_mesh;
  SP_state d_state;
  SP_quadrature d_quadrature;
  SP_externalsource d_externalsource;
  SP_fissionsource d_fissionsource;
  std::string d_problem_type;

  /// \}

  /// \name Implementation
  /// \{

  /// Generic problem setup.
  void setup();

  /// \}

};

template <class D>
PyExecute<D>::PyExecute(SP_input input,
                        SP_material material,
                        SP_mesh mesh)
 : d_input(input)
 , d_material(material)
 , d_mesh(mesh)
{
  // Check preconditions.
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Setup the problem.
  setup();

  // Check postconditions.
  Ensure(d_quadrature);
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

  //-------------------------------------------------------------------------//
  // Quadrature
  //-------------------------------------------------------------------------//

  string quad_type;
  if (!d_input->check("quad_type"))
  {
    if (D::dimension == 1) quad_type = "gausslegendre";
    if (D::dimension == 2) quad_type = "quadruplerange";
    if (D::dimension == 3) quad_type = "levelsymmetric";
  }
  else
  {
    quad_type = d_input->get<string>("quad_type");
  }
  int quad_order;
  if (!d_input->check("quad_order"))
  {
    quad_order = 2;
  }
  else
  {
    quad_order = d_input->get<int>("quad_order");
  }
  QuadratureFactory quad_factory;
  quad_factory.build(d_quadrature, quad_type, quad_order, D::dimension);

  //-------------------------------------------------------------------------//
  // State
  //-------------------------------------------------------------------------//

  d_state = new State(d_input, d_mesh, d_quadrature);

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

  typename Boundary<D>::SP_boundary boundary;
  boundary = new Boundary<D>(d_input, d_mesh, d_quadrature);

  //--------------------------------------------------------------------------//
  // Create solver and solve.
  //--------------------------------------------------------------------------//

  if (d_problem_type == "eigenvalue")
  {
    Require(d_fissionsource);
    PowerIteration<D> solver(d_input,
                             d_state,
                             d_mesh,
                             d_material,
                             d_quadrature,
                             boundary,
                             d_externalsource,
                             d_fissionsource);
    solver.solve();
  }
  else if (d_problem_type == "fixed" || d_problem_type == "fixed_multiply")
  {
    GaussSeidel<D> solver(d_input,
                          d_state,
                          d_mesh,
                          d_material,
                          d_quadrature,
                          boundary,
                          d_externalsource,
                          d_fissionsource);
    solver.solve();
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
