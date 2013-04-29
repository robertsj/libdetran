//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Execute.cc
 *  @author robertsj
 *  @date   Apr 11, 2012
 *  @brief  Execute class definition.
 */
//---------------------------------------------------------------------------//

#include "Execute.hh"
#include "external_source/ConstantSource.hh"
#include "external_source/DiscreteSource.hh"
#include "external_source/IsotropicSource.hh"
#include "utilities/Definitions.hh"
#ifdef DETRAN_ENABLE_SILO
#include "ioutils/SiloOutput.hh"
#endif

namespace detran
{

Execute::Execute(StupidParser &parser)
 : d_dimension(0)
{
  // Parse the input.
  d_input = parser.parse_input();
  // Parse the material.
  d_material = parser.parse_material();
  // Parse the mesh.
  d_mesh = parser.parse_mesh();
  Assert(d_input);
  Assert(d_material);
  Assert(d_mesh);

  // Get the dimension of the problem.
  d_dimension = d_input->get<int>("dimension");
  // Problem type
  if (d_input->check("problem_type"))
    d_problem_type = d_input->get<std::string>("problem_type");
  // Number of groups
  d_number_groups = d_input->get<int>("number_groups");

  // Postconditions
  Ensure(d_dimension > 0 && d_dimension < 4);
}

template<class D>
void Execute::solve()
{

  using std::string;
  using detran_utilities::vec_dbl;
  using detran_utilities::vec_int;
  typename FixedSourceManager<D>::SP_quadrature q;

  //--------------------------------------------------------------------------//
  // FIXED SOURCE
  //--------------------------------------------------------------------------//

  if (d_problem_type != "eigenvalue")
  {

    // Must have a RHS for fixed source problems
    if (!d_input->check("source_type"))
    {
      std::cout << "No source type is listed in the input; skipping solve."
                << std::endl;
    }

    // Set multiplying fixed source problem flag and get the scaling factor
    bool multiply = false;
    if (d_problem_type == "multiply") multiply = true;
    double keff = 1.0;
    if (d_input->check("fixed_scaling_factor"))
      keff = d_input->get<double>("fixed_scaling_factor");

    // Create fixed source solver and setup
    FixedSourceManager<D> solver(d_input, d_material, d_mesh, multiply);
    solver.setup();

    // Get the quadrature (may be NULL)
    q = solver.quadrature();

    // Define the fixed source
    string source_type = d_input->get<string>("source_type");
    if (source_type == "constant")
    {

      double strength = 1.0;
      if (d_input->check("source_type"))
        strength = d_input->get<double>("source_strength");
      d_externalsource = new detran_external_source::
        ConstantSource(d_number_groups, d_mesh, strength, q);

    }
    else if (source_type == "isotropic")
    {

      Insist(d_input->check("source_spectra"),
             "Isotropic source requested, but no source spectra provided!");
      vec_dbl spectrav = d_input->get<vec_dbl>("source_spectra");
      int number_spectra = spectrav.size() / d_number_groups;
      detran_external_source::IsotropicSource::spectra_type
        spectra(number_spectra, vec_dbl(d_number_groups));
      int count = 0;
      for (int n = 0; n < number_spectra; n++)
        for (int g = 0; g < d_number_groups; g++)
          spectra[n][g] = spectrav[count++];
      Insist(d_input->check("source_spectra_map"),
             "Isotropic source requested, but no source map provided!");
      vec_int spectra_map = d_input->get<vec_int>("source_spectra_map");
      d_externalsource = new detran_external_source::
        IsotropicSource(d_number_groups, d_mesh, spectra, spectra_map, q);

    }
    else
    {
      THROW("Unsupported external source requested: " + source_type);
    }

    // Add the source.  Limited to one source via this interface.
    solver.set_source(d_externalsource);

    // Solve the problem
    solver.solve(keff);

    // Get the state and quadrature
    d_state = solver.state();
    q = solver.quadrature();
  }

  //--------------------------------------------------------------------------//
  // EIGENVALUE
  //--------------------------------------------------------------------------//

  else if (d_problem_type == "eigenvalue")
  {
    // Create eigenvalue solver source solver
    EigenvalueManager<D> solver(d_input, d_material, d_mesh);

    // Solve the problem
    solver.solve();

    // Get the state
    d_state = solver.state();

    // Spit out the eigenvalue for sanity (need better output)
    std::cout << " *** Eigenvalue = " << d_state->eigenvalue() << std::endl;
  }

  else
  {
    std::cout << "Unsupported problem type given.  Options are:" << std::endl
              << " -- fixed" << std::endl
              << " -- multiply" << std::endl
              << " -- eigenvalue" << std::endl
              << std::endl;
  }

  //--------------------------------------------------------------------------//
  // OUTPUT
  //--------------------------------------------------------------------------//

  if (d_input->check("output_silo"))
  {

#ifdef DETRAN_ENABLE_SILO
    // Create the silo output
    detran_ioutils::SiloOutput silo(d_mesh);
    std::string siloname = d_input->get<std::string>("output_silo");
    silo.initialize(siloname);

    // Write out mesh maps; some are optional.
    silo.write_mesh_map("MATERIAL");
    silo.write_mesh_map("PINS");
    silo.write_mesh_map("ASSEMBLIES");

    // Write scalar flux
    silo.write_scalar_flux(d_state);

    // Write angular flux
    if (d_state->store_angular_flux())
    {
      silo.write_angular_flux(d_state, q);
    }

    silo.finalize();
#else
    std::cout << "Silo output requested, but Silo is unavaible." << std::endl;
#endif
  }

}

// Explicit instantiations of function templates 
template void Execute::solve<_1D>();
template void Execute::solve<_2D>();
template void Execute::solve<_3D>();

} // end namespace detran


