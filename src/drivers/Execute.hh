//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Execute.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  Execute class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EXECUTE_HH_
#define EXECUTE_HH_

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
#include <string>

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

  /// Do nontemplated portion of problem setup.
  void setup();

};

template<class D>
void Execute::solve()
{

  //--------------------------------------------------------------------------//
  // Boundary
  //--------------------------------------------------------------------------//

  typename Boundary<D>::SP_boundary boundary;
  Assert(d_quadrature);
  boundary = new Boundary<D>(d_input, d_mesh, d_quadrature);

  //--------------------------------------------------------------------------//
  // Create solver and solve.
  //--------------------------------------------------------------------------//

  // \todo Perhaps for outers, use a setup with a generic interface.

  if (d_problem_type == "eigenvalue")
  {
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
    THROW("Unsupported problem type given.");
  }

  State::moments_type phi = d_state->phi(0);

  std::cout << phi[0] << " " << phi[1] << std::endl;

}

} // end namespace detran

#endif /* EXECUTE_HH_ */

//---------------------------------------------------------------------------//
//              end of Execute.hh
//---------------------------------------------------------------------------//
