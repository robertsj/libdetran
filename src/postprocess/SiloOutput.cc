//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SiloOutput.cc
 * \brief  SiloOutput 
 * \author Jeremy Roberts
 * \date   Jul 27, 2012
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_SILO

// Postprocess
#include "SiloOutput.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"

// System
#include <string>
#include <stdio.h>

namespace detran_postprocess
{

SiloOutput::SiloOutput(SP_input input, SP_mesh mesh)
  : d_input(input)
  , d_mesh(mesh)
  , d_initialized(false)
{
  // Preconditions
  Require(d_input);
  Require(d_mesh);
}

SiloOutput::~SiloOutput()
{
  if(d_initialized) finalize();
}

void SiloOutput::initialize()
{
  if (!d_initialized)
  {

    // Get silo name from input.
    std::string filename = "detran.silo";
    if (d_input->check("silo_filename"))
    {
      filename = d_input->get<std::string>("silo_filename");
    }

    // Create the file
    d_silofile = DBCreate((char *)filename.c_str(), DB_CLOBBER,
                          DB_LOCAL, NULL, DB_PDB );
    d_initialized = true;

    // Check for errors
    Insist(d_silofile, "Error creating silo file.");
  }
}

void SiloOutput::finalize()
{
  if (d_initialized)
  {
    DBClose(d_silofile);
    d_initialized = false;
  }
}

void SiloOutput::write_flux(SP_state state)
{
  using std::cout;
  using std::endl;

  if (!d_initialized)
  {
    std::cout << "Not initialized; not writing fluxes" << std::endl;
    return;
  }

  // Problem dimension
  int dimension = d_mesh->dimension();
  if (dimension == 1) dimension = 2;

  // Names of the coordinates
  char *coordnames[dimension];

  // Array of coordinate arrays
  double *coords[dimension];

  // Number of nodes in each dimension
  int dims[dimension];

  // Name the coordinates
  char c1[5] = "x";
  coordnames[0] = c1;
  char c2[5] = "y";
  coordnames[1] = c2;
  if (dimension == 3)
  {
    char c3[5] = "z";
    coordnames[1] = c3;
  }

  // Set coordinate dimensions
  dims[0] = d_mesh->number_cells_x() + 1;
  dims[1] = d_mesh->number_cells_y() + 1;
  if (dimension == 3) dims[2] = d_mesh->number_cells_z() + 1;

  // Create the coordinate vectors (mesh edges)
  detran::vec_dbl x(d_mesh->number_cells_x() + 1, 0.0);
  detran::vec_dbl y(d_mesh->number_cells_y() + 1, 0.0);
  detran::vec_dbl z(d_mesh->number_cells_z() + 1, 0.0);
  for (int i = 0; i < d_mesh->number_cells_x(); i++)
  {
    x[i + 1] = x[i] + d_mesh->dx(i);
  }
  for (int i = 0; i < d_mesh->number_cells_y(); i++)
    y[i + 1] = y[i] + d_mesh->dy(i);
  for (int i = 0; i < d_mesh->number_cells_z(); i++)
    z[i + 1] = z[i] + d_mesh->dz(i);
  coords[0] = &x[0];
  coords[1] = &y[0];
  if (dimension == 3) coords[2] = &z[0];

  // Add the mesh
  DBPutQuadmesh(d_silofile, "mesh", NULL, coords, dims, 2, DB_DOUBLE, DB_COLLINEAR, NULL);

  for (int d = 0; d < dimension; d++)
    dims[d]--;

  // Get the scalar values from the fluxes
  Assert(d_input->check("number_groups"));
  for (int g = 0; g < d_input->get<int>("number_groups"); g++)
  {

    // Write the group
    char buffer[14];
    sprintf(buffer, "group%i", g);

    // Get a pointer to the group flux
    double *phi = &state->phi(g)[0];

    // Write to silo
    DBPutQuadvar1(d_silofile, buffer, "mesh", phi,
                  dims, dimension, NULL, 0, DB_DOUBLE,
                  DB_ZONECENT, NULL);

  }
}

} // end namespace detran_postprocess


#endif // DETRAN_ENABLE_SILO

//---------------------------------------------------------------------------//
//              end of file SiloOutput.cc
//---------------------------------------------------------------------------//
