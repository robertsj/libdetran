//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Execute.cc
 * \author robertsj
 * \date   Apr 11, 2012
 * \brief  Execute class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "Execute.hh"

// Quadratures.
#include "QuadrupleRange.hh"
#include "GaussLegendre.hh"
#include "QuadratureFactory.hh"

// MOC
#include "Tracker.hh"

// Sources
#include "ConstantSource.hh"
#include "IsotropicSource.hh"

// Solvers
#include "PowerIteration.hh"
#include "GaussSeidel.hh"

// Utilities
#include "Definitions.hh"

namespace detran
{

Execute::Execute(StupidParser &parser)
 : d_dimension(0)
{
  // Parse the input.
  d_input = parser.parse_input();
  Ensure(d_input);

  // Parse the material.
  d_material = parser.parse_material();
  Ensure(d_material);

  // Parse the mesh.
  d_mesh = parser.parse_mesh();
  Ensure(d_mesh);

  // Setup other things.
  setup();

  // Postconditions
  Ensure(d_dimension);
}

//--------------------------------------------------------------------------//
// Private implementation.
//--------------------------------------------------------------------------//

void Execute::setup()
{
  using std::string;

  // Get the dimension of the problem.
  d_dimension = d_input->get<int>("dimension");

  // Number of groups
  int number_groups = d_input->get<int>("number_groups");

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

  //--------------------------------------------------------------------------//
  // Quadrature
  //--------------------------------------------------------------------------//

  QuadratureFactory quad_factory;
  quad_factory.build(d_quadrature, d_input, d_dimension);
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

  //--------------------------------------------------------------------------//
  // External source
  //--------------------------------------------------------------------------//

  if (d_input->check("source_type") and d_problem_type != "eigenvalue")
  {
    string source_type = d_input->get<string>("source_type");
    if (source_type == "constant")
    {
      double strength = 1.0;
      if (d_input->check("source_type"))
        strength = d_input->get<double>("source_strength");
      d_externalsource =
          new ConstantSource(d_mesh, d_quadrature, number_groups, strength);
    }
    else if (source_type == "isotropic")
    {
      Insist(d_input->check("source_spectra"),
             "Isotropic source requested, but no source spectra provided!");
      vec_dbl spectrav = d_input->get<vec_dbl>("source_spectra");
      int number_spectra = spectrav.size() / number_groups;
      IsotropicSource::spectra_type spectra(number_spectra, vec_dbl(number_groups));
      int count = 0;
      for (int n = 0; n < number_spectra; n++)
        for (int g = 0; g < number_groups; g++)
          spectra[n][g] = spectrav[count++];
      Insist(d_input->check("source_spectra_map"),
             "Isotropic source requested, but no source map provided!");
      vec_int spectra_map =  d_input->get<vec_int>("source_spectra_map");
      d_externalsource = new IsotropicSource(d_mesh, d_quadrature,
                                             number_groups, spectra,
                                             spectra_map);
    }
    else
    {
      THROW("Unsupported external source requested: " + source_type);
    }
  }

  //--------------------------------------------------------------------------//
  // State
  //--------------------------------------------------------------------------//

  d_state = new State(d_input, d_mesh, d_quadrature);

  //--------------------------------------------------------------------------//
  // Fission source
  //--------------------------------------------------------------------------//

  if (d_problem_type == "eigenvalue" || d_problem_type == "fixed_multiply")
  {
    d_fissionsource = new FissionSource(d_state, d_mesh, d_material);
  }

}


} // end namespace detran


