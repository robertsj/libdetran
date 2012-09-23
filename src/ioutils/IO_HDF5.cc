//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   IO_HDF5.cc
 * \brief  IO_HDF5 
 * \author Jeremy Roberts
 * \date   Jul 29, 2012
 * \todo   There is a lot of error checking that could be done more smoothly
 * \todo   There is a lot of refactoring to do to reduce bulk
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.hh"

#ifdef DETRAN_ENABLE_HDF5

#include "IO_HDF5.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include <sstream>

namespace detran_ioutils
{

IO_HDF5::IO_HDF5(std::string filename)
  : d_filename(filename)
  , d_open(false)
{
  /* ... */
}

void IO_HDF5::open()
{
  // Open the HDF5 file
  d_file_id = H5Fcreate(d_filename.c_str(), // filename
                        H5F_ACC_TRUNC,      // overwrite existing file
                        H5P_DEFAULT,        // file create property list
                        H5P_DEFAULT);       // file access property list
  d_open = true;
}

void IO_HDF5::write(SP_input input)
{
  // Preconditions
  Require(input);

  if (!d_open) open();

  // Create the input group
  hid_t group = H5Gcreate(d_file_id, "input",
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write the data.
  write_map(group, "int_data", input->get_map<int>());
  write_map(group, "dbl_data", input->get_map<double>());
  write_map(group, "str_data", input->get_map<std::string>());
  write_map(group, "vec_int_data", input->get_map<vec_int>());
  write_map(group, "vec_dbl_data", input->get_map<vec_dbl>());

  // Close the group.
  herr_t status = H5Gclose(group);
}

void IO_HDF5::write(SP_material mat)
{
  // Preconditions
  Require(mat);

  if (!d_open) open();

  // Create the material group
  hid_t group = H5Gcreate(d_file_id, "material",
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status;

  //-------------------------------------------------------------------------//
  // ATTRIBUTES
  //-------------------------------------------------------------------------//

  hid_t att_space;
  hid_t att;

  // Create scalar attribute.
  att_space = H5Screate(H5S_SCALAR);
  att = H5Acreate(group, "number_groups", H5T_NATIVE_INT, att_space,
                  H5P_DEFAULT, H5P_DEFAULT);
  // Write scalar attribute.
  int ng = mat->number_groups();
  status = H5Awrite(att, H5T_NATIVE_INT, (void *) &ng);
  // Close dataspace and attribute
  status = H5Sclose(att_space);
  status = H5Aclose(att);

  // Create scalar attribute.
  att_space = H5Screate(H5S_SCALAR);
  att = H5Acreate(group, "number_materials", H5T_NATIVE_INT, att_space,
                  H5P_DEFAULT, H5P_DEFAULT);
  // Write scalar attribute.
  int nm = mat->number_materials();
  status = H5Awrite(att, H5T_NATIVE_INT, (void *) &nm);
  // Close dataspace and attribute
  status = H5Sclose(att_space);
  status = H5Aclose(att);

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  hid_t dset;
  hid_t space;

  for (int m = 0; m < nm; m++)
  {
    // Create material name
    std::string name = "material";
    std::ostringstream convert;
    convert << m;
    name += convert.str();

    // Create the material group
    hid_t group_m = H5Gcreate(group, name.c_str(),
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    // All data is ng long, except scatter, which is ng * ng
    hsize_t dims1[1] = {ng};
    hsize_t dims2[2] = {ng, ng};

    // TOTAL, FISSION, NU, CHI, DIFFUSION
    space  = H5Screate_simple(1, dims1, NULL);

    dset   = H5Dcreate(group_m, "sigma_t", H5T_NATIVE_DOUBLE, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &mat->sigma_t(m)[0]);
    status = H5Dclose(dset);

    dset   = H5Dcreate(group_m, "sigma_f", H5T_NATIVE_DOUBLE, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &mat->sigma_f(m)[0]);
    status = H5Dclose(dset);

    dset   = H5Dcreate(group_m, "nu", H5T_NATIVE_DOUBLE, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &mat->nu(m)[0]);
    status = H5Dclose(dset);

    dset   = H5Dcreate(group_m, "chi", H5T_NATIVE_DOUBLE, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &mat->chi(m)[0]);
    status = H5Dclose(dset);

    dset   = H5Dcreate(group_m, "diff_coef", H5T_NATIVE_DOUBLE, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &mat->diff_coef(m)[0]);
    status = H5Dclose(dset);

    // SCATTER
    space  = H5Screate_simple(2, dims2, NULL);

    // There may be a better way to avoid so much copying.  Even so, this
    // is *not* intended to be an interface for large data.
    double *ss;
    ss = new double[ng * ng];
    for (int g = 0; g < ng; g++)
      for (int gp = 0; gp < ng; gp++)
        ss[gp + g * ng] = mat->sigma_s(m, g, gp);

    dset   = H5Dcreate(group_m, "sigma_s", H5T_NATIVE_DOUBLE, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, ss);
    status = H5Dclose(dset);
    delete [] ss;

    status = H5Sclose(space);

    // Close the group.
    status = H5Gclose(group_m);
  }

  // Close the group.
  status = H5Gclose(group);

}

void IO_HDF5::write(SP_mesh mesh)
{
  // Preconditions
  Require(mesh);

  if (!d_open) open();

  // Create the mesh group
  hid_t group = H5Gcreate(d_file_id, "mesh",
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write the required attributes.
  write_scalar_attribute(group, "dimension", mesh->dimension());
  write_scalar_attribute(group, "number_cells_x", mesh->number_cells_x());
  write_scalar_attribute(group, "number_cells_y", mesh->number_cells_y());
  write_scalar_attribute(group, "number_cells_z", mesh->number_cells_z());

  // Write the dx, dy, and dz vectors.
  bool status;
  status = write_vec(group, "dx", mesh->dx());
  Assert(status);
  status = write_vec(group, "dy", mesh->dy());
  Assert(status);
  status = write_vec(group, "dz", mesh->dz());
  Assert(status);

  // Write all mesh maps.  Note, at *least* MATERIAL must be present.
  detran_geometry::Mesh::mesh_map_type map = mesh->get_mesh_map();

  status = write_map(group, "mesh_map", map);
  Assert(status);
}

void IO_HDF5::close()
{
  // If not open, ignore.
  if (!d_open) return;

  // Close the HDF5 file
  H5Fclose(d_file_id);
  d_open = false;
}

detran_utilities::InputDB::SP_input IO_HDF5::read_input()
{
  // Preconditions
  /* ... */

  // Create the input object.
  SP_input input(new detran_utilities::InputDB());

  // Open the file if necessary.
  if (!d_open)
    d_file_id = H5Fopen(d_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the input group
  htri_t flag = H5Lexists(d_file_id, "input", H5P_DEFAULT);
  Insist(flag > 0, "Input group does not exist; can't read data.");
  hid_t group = H5Gopen(d_file_id, "input", H5P_DEFAULT);

  // Read data
  read_data<int>(input, group, "int_data");
  read_data<double>(input, group, "dbl_data");
  read_data<std::string>(input, group, "str_data");
  read_data<vec_int>(input, group, "vec_int_data");
  read_data<vec_dbl>(input, group, "vec_dbl_data");

//  read_map(group, "int_data", input->get_map<int>());
//  read_map(group, "dbl_data", input->get_map<double>());
//  read_map(group, "str_data", input->get_map<std::string>());
//  read_map(group, "vec_int_data", input->get_map<vec_int>());
//  read_map(group, "vec_dbl_data", input->get_map<vec_dbl>());

  // Close the group.
  herr_t status = H5Gclose(group);

  // Postconditions
  Ensure(input);
  return input;
}

detran_material::Material::SP_material IO_HDF5::read_material()
{
  // Preconditions
  /* ... */

  // Material to be filled.
  SP_material mat;

  // Open the file if necessary.
  if (!d_open)
    d_file_id = H5Fopen(d_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // HDF5 bool
  htri_t flag;
  herr_t status;

  // Check if the material group exists.
  flag = H5Lexists(d_file_id, "material", H5P_DEFAULT);
  Insist(flag > 0, "Material group does not exist; can't read data.");
  hid_t group = H5Gopen(d_file_id, "material", H5P_DEFAULT);

  //-------------------------------------------------------------------------//
  // ATTRIBUTES
  //-------------------------------------------------------------------------//

  // Read them.
  int ng;
  int nm;

  Insist(read_scalar_attribute(group, "number_groups", ng),
    "Number of groups missing from HDF5 file.");

  Insist(read_scalar_attribute(group, "number_materials", nm),
    "Number of materials missing from HDF5 file.");

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Create the material object.
  mat = new detran_material::Material(nm, ng, false);

  hid_t dset;
  hid_t space;

  // All data is ng long, except scatter, which is ng * ng
  hsize_t dims1[1] = {ng};
  hsize_t dims2[2] = {ng, ng};

  for (int m = 0; m < nm; m++)
  {
    // Create material name
    std::string name = "material";
    std::ostringstream convert;
    convert << m;
    name += convert.str();

    // Switch to this material's group
    hid_t group_m = H5Gopen(group, name.c_str(), H5P_DEFAULT);


    // TOTAL, FISSION, NU, CHI, DIFFUSION

    // Read buffer.
    vec_dbl v(ng, 0.0);

    Insist(read_vec(group_m, "sigma_t", v),
      "Error reading SigmaT from HDF5.  SigmaT is *required*");
    mat->set_sigma_t(m, v);

    if (read_vec(group_m, "sigma_f", v))
      mat->set_sigma_f(m, v);
    if (read_vec(group_m, "nu", v))
      mat->set_nu(m, v);
    if (read_vec(group_m, "chi", v))
      mat->set_chi(m, v);
    if (read_vec(group_m, "diff_coef", v))
      mat->set_diff_coef(m, v);

    // SCATTER

    // Read buffer
    double *ss;
    ss = new double[ng * ng];
    dset   = H5Dopen(group_m, "sigma_s", H5P_DEFAULT);
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, ss);
    status = H5Dclose(dset);

    for (int g = 0; g < ng; g++)
      for (int gp = 0; gp < ng; gp++)
        mat->set_sigma_s(m, g, gp, ss[gp + g * ng]);
    delete [] ss;

    // Close the this material's group.
    status = H5Gclose(group_m);
  }

  mat->finalize();

  // Close the group.
  status = H5Gclose(group);

  // Postconditions
  /* ... */

  return mat;
}

detran_geometry::Mesh::SP_mesh IO_HDF5::read_mesh()
{
  // Preconditions
  /* ... */

  SP_mesh mesh;

  // Open the file if necessary.
  if (!d_open)
    d_file_id = H5Fopen(d_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  d_open = true;

  // Open the input group
  Insist(exists(d_file_id, "mesh") > 0,
    "Mesh group does not exist; can't read data.");
  hid_t group = H5Gopen(d_file_id, "mesh", H5P_DEFAULT);

  // Read required attributes.
  int dim, nx, ny, nz;
  Insist(read_scalar_attribute(group, "dimension", dim),
    "Mesh dimension attribute missing from HDF5 file.");
  Insist(read_scalar_attribute(group, "number_cells_x", nx),
    "Number of x cells attribute missing from HDF5 file.");
  Insist(read_scalar_attribute(group, "number_cells_y", ny),
    "Number of y cells attribute missing from HDF5 file.");
  Insist(read_scalar_attribute(group, "number_cells_z", nz),
    "Number of z cells attribute missing from HDF5 file.");
  Insist(nx > 0 and ny > 0 and nz > 0,
    "The number of cells per direction must be positive");
  Insist(dim == 1 or dim == 2 or dim == 3,
    "The mesh dimension attribute must be 1, 2, or 3");

  // Read the discretization data.
  vec_dbl dx(nx, 0.0);
  vec_dbl ex(nx + 1, 0.0);
  Insist(read_vec(group, "dx", dx), "Failed to read dx from HDF5 file.");
  for (int i = 0; i < nx; i++)
    ex[i + 1] = ex[i] + dx[i];
  vec_dbl dy(ny, 0.0);
  vec_dbl ey(ny + 1, 0.0);
  Insist(read_vec(group, "dy", dy), "Failed to read dx from HDF5 file.");
  for (int i = 0; i < ny; i++)
    ey[i + 1] = ey[i] + dy[i];
  vec_dbl dz(nz, 0.0);
  vec_dbl ez(nz + 1, 0.0);
  Insist(read_vec(group, "dz", dz), "Failed to read dx from HDF5 file.");
  for (int i = 0; i < nz; i++)
    ez[i + 1] = ez[i] + dz[i];

  // Fill the mesh maps, and extract the material map.
  detran_geometry::Mesh::mesh_map_type map;
  Insist(read_map(group, "mesh_map", map),
    "Problem reading mesh map.  It is missing or empty.");
  vec_int mt = map["MATERIAL"];

  // Create the mesh.
  if (dim == 1)
    mesh = new detran_geometry::Mesh1D(ex, mt);
  else if (dim == 2)
    mesh = new detran_geometry::Mesh2D(ex, ey, mt);
  else
    mesh = new detran_geometry::Mesh3D(ex, ey, ez, mt);

  // Add the rest of the maps.
  detran_geometry::Mesh::mesh_map_type::iterator it = map.begin();
  for (; it != map.end(); it++)
  {
    if (it->first != "MATERIAL")
      mesh->add_mesh_map(it->first, it->second);
  }

  // Close the group.
  herr_t status = H5Gclose(group);

  // Postconditions
  Ensure(mesh);

  return mesh;
}


} // end namespace detran_ioutils

#endif // DETRAN_ENABLE_HDF5

//---------------------------------------------------------------------------//
//              end of file IO_HDF5.cc
//---------------------------------------------------------------------------//
