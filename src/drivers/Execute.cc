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

  Require(d_quadrature);
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

  //--------------------------------------------------------------------------//
  // Quadrature
  //--------------------------------------------------------------------------//

  string quad_type;
  if (!d_input->check("quad_type"))
  {
    if (d_dimension == 1) quad_type = "gausslegendre";
    if (d_dimension == 2) quad_type = "quadruplerange";
    if (d_dimension == 3) quad_type = "levelsymmetric";
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
  quad_factory.build(d_quadrature, quad_type, quad_order, d_dimension);
  Require(d_quadrature);

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
      THROW("Unsupported external source requested.")
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


