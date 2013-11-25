//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   IO_HDF5.cc
 *  @brief  IO_HDF5 member definitions
 *  @author Jeremy Roberts
 *  @date   Jul 29, 2012
 *  @todo   There is a lot of error checking that could be done more smoothly
 *  @todo   There is a lot of refactoring to do to reduce bulk
 */
//---------------------------------------------------------------------------//

#include "detran_config.hh"
#include "IO_HDF5.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include <sstream>

namespace detran_ioutils
{

#ifdef DETRAN_ENABLE_HDF5

//---------------------------------------------------------------------------//
IO_HDF5::IO_HDF5(std::string filename)
  : d_filename(filename)
  , d_open(false)
  , d_file_id(0)
{
  /* ... */
}

//---------------------------------------------------------------------------//
IO_HDF5::~IO_HDF5()
{
  if (d_open) close();
}

//---------------------------------------------------------------------------//
void IO_HDF5::open(const int flag)
{
  Require(flag < END_HDF5_FILE_ACCESS);
  const char* fname = d_filename.c_str();
  if (flag == HDF5_READ_ONLY)
    d_file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  else if (flag == HDF5_OVERWRITE)
    d_file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  Insist(d_file_id > 0, "Error opening HDF5 file " + d_filename);
  d_open = true;
}

//---------------------------------------------------------------------------//
void IO_HDF5::write(SP_input input)
{
  Require(input);
  if (!d_open) open(HDF5_OVERWRITE);
  write(input, "input", d_file_id);
}

//---------------------------------------------------------------------------//
void IO_HDF5::write(SP_input input, std::string name, hid_t root)
{
  Require(input);
  Require(d_open);

  // Create the input group
  hid_t group = H5Gcreate(root, name.c_str(),
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write the data.
  write_map(group, "int_data",     input->get_map<int>());
  write_map(group, "dbl_data",     input->get_map<double>());
  write_map(group, "str_data",     input->get_map<std::string>());
  write_map(group, "vec_int_data", input->get_map<vec_int>());
  write_map(group, "vec_dbl_data", input->get_map<vec_dbl>());
  write_map(group, "db_data",      input->get_map<SP_input>());

  // Close the group.
  herr_t status = H5Gclose(group);
}

//---------------------------------------------------------------------------//
void IO_HDF5::write(SP_material mat)
{
  Require(mat);

  if (!d_open) open(HDF5_OVERWRITE);

  // Create the material group
  hid_t group = H5Gcreate(d_file_id, "material",
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  //-------------------------------------------------------------------------//
  // ATTRIBUTES
  //-------------------------------------------------------------------------//

  bool flag;

  int ng = mat->number_groups();
  flag = write_scalar_attribute(group, "number_groups", ng);
  Assert(flag);

  int nm = mat->number_materials();
  flag = write_scalar_attribute(group, "number_materials", nm);
  Assert(flag);

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  herr_t status;
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
    hsize_t dims2[2] = {ng, ng};

    // TOTAL, FISSION, NU, CHI, DIFFUSION
    write_vec(group_m, "sigma_t",     mat->sigma_t(m));
    write_vec(group_m, "sigma_a",     mat->sigma_a(m));
    write_vec(group_m, "sigma_f",     mat->sigma_f(m));
    write_vec(group_m, "nu",          mat->nu(m));
    write_vec(group_m, "chi",         mat->chi(m));
    write_vec(group_m, "diff_coef",   mat->sigma_a(m));

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

    status = H5Gclose(group_m);
  }

  // Close the group.
  status = H5Gclose(group);
}

//---------------------------------------------------------------------------//
void IO_HDF5::write(SP_mesh mesh)
{
  // Preconditions
  Require(mesh);

  if (!d_open) open(HDF5_OVERWRITE);

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

//---------------------------------------------------------------------------//
void IO_HDF5::close()
{
  // If not open, ignore.
  if (!d_open) return;

  // Close the HDF5 file
  H5Fclose(d_file_id);
  d_open = false;
}

//---------------------------------------------------------------------------//
IO_HDF5::SP_input IO_HDF5::read_input()
{
  // Open the file if necessary.
  if (!d_open) open(HDF5_READ_ONLY);

  // Create the database
  SP_input db = read_input(d_file_id, "input");

  // Postconditions
  Ensure(db);
  return db;
}

//---------------------------------------------------------------------------//
IO_HDF5::SP_input IO_HDF5::read_input(hid_t root, const char* name)
{
  Require(d_open);

  // Create the input object.
  std::string str_name(name);
  SP_input db(new detran_utilities::InputDB(str_name));

  // Open the input group
  htri_t flag = H5Lexists(root, name, H5P_DEFAULT);
  Insist(flag > 0, "Group " + str_name + " does not exist; can't read data.");
  hid_t group = H5Gopen(root, name, H5P_DEFAULT);

  // Read data
  read_data<int>(db, group,           "int_data");
  read_data<double>(db, group,        "dbl_data");
  read_data<std::string>(db, group,   "str_data");
  read_data<vec_int>(db, group,       "vec_int_data");
  read_data<vec_dbl>(db, group,       "vec_dbl_data");
  read_data<SP_input>(db, group,      "db_data");

  // Close the group.
  herr_t status = H5Gclose(group);

  Ensure(db);
  return db;
}

//---------------------------------------------------------------------------//
detran_material::Material::SP_material IO_HDF5::read_material()
{
  // Material to be filled.
  SP_material mat;

  // Open the file if necessary.
  if (!d_open) open(HDF5_READ_ONLY);

  // HDF5 bool
  htri_t flag;
  herr_t status;

  // Check if the material group exists.
  Insist(exists(d_file_id, "material"),
    "Material group does not exist; can't read data.");
  hid_t group = H5Gopen(d_file_id, "material", H5P_DEFAULT);

  //-------------------------------------------------------------------------//
  // ATTRIBUTES
  //-------------------------------------------------------------------------//

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
  mat = new detran_material::Material(nm, ng, d_filename + "_material");

  hid_t dset;
  hid_t space;

  for (int m = 0; m < nm; m++)
  {
    // Create material name
    std::string name = "material";
    std::ostringstream convert;
    convert << m;
    name += convert.str();

    // Switch to this material's group
    Insist(exists(group, name.c_str()), "Group not found");
    hid_t group_m = H5Gopen(group, name.c_str(), H5P_DEFAULT);

    // TOTAL, FISSION, NU, CHI, DIFFUSION
    vec_dbl v(ng, 0.0);
    Insist(read_vec(group_m, "sigma_t", v),
      "Error reading SigmaT from HDF5.  SigmaT is *required*");
    mat->set_sigma_t(m, v);
    if (read_vec(group_m, "sigma_a", v))    mat->set_sigma_a(m, v);
    if (read_vec(group_m, "sigma_f", v))    mat->set_sigma_f(m, v);
    if (read_vec(group_m, "nu", v))         mat->set_nu(m, v);
    if (read_vec(group_m, "chi", v))        mat->set_chi(m, v);
    if (read_vec(group_m, "diff_coef", v))  mat->set_diff_coef(m, v);

    // SCATTER
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

  return mat;
}

//---------------------------------------------------------------------------//
detran_geometry::Mesh::SP_mesh IO_HDF5::read_mesh()
{
  // Mesh to fill
  SP_mesh mesh;

  // Open the file if necessary.
  if (!d_open) open(HDF5_READ_ONLY);

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
  Insist(nx > 0 && ny > 0 && nz > 0,
    "The number of cells per direction must be positive");
  Insist(dim == 1 || dim == 2 || dim == 3,
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

  Ensure(mesh);
  return mesh;
}

//---------------------------------------------------------------------------//
herr_t find_db_groups(hid_t              loc_id,
                      const char        *name,
                      const H5L_info_t  *info,
                      void              *data)
{
  herr_t status;
  H5O_info_t infobuf;
  std::vector<std::string> &names = *((std::vector<std::string>*)data);

  // Get type of the object.  There should *only* be groups.
  status = H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
  Insist(infobuf.type == H5O_TYPE_GROUP,
         "Something other than a group was found nested in a parameter DB.");

  // Keep the name
  names.push_back(std::string(name));

  return 0;
}

//---------------------------------------------------------------------------//
template <>
bool IO_HDF5::read_data<IO_HDF5::SP_input>(SP_input      db,
                                           hid_t         root,
                                           std::string   name)
{
  Require(db);

  // Open the root db group
  if (!exists(root, name.c_str())) return false;
  hid_t group = H5Gopen(root, "db_data", H5P_DEFAULT);

  // Note, we need to iterate through the groups in ~/db_data, as these
  // represent imbedded db's.
  std::vector<std::string> names;
  hid_t status = H5Literate(group,
                            H5_INDEX_NAME,
                            H5_ITER_NATIVE,
                            NULL,
                            find_db_groups,
                            (void *) &names);

  // Loop through the nested groups, and add the resulting db's to the main db.
  for (int i = 0; i < names.size(); ++i)
  {
    Assert(exists(group, names[i].c_str()));
    SP_input nested_db = read_input(group, names[i].c_str());
    db->put<SP_input>(names[i], nested_db);
  }

  return true;
}

#else

//-------------------------------------------------------------------------//
IO_HDF5::IO_HDF5(std::string filename)
{
  THROW("NOT IMPLEMENTED");
}

//-------------------------------------------------------------------------//
IO_HDF5::~IO_HDF5(){}

//-------------------------------------------------------------------------//
void IO_HDF5::open(const int flag){}

//-------------------------------------------------------------------------//
void IO_HDF5::write(SP_input input){}

//-------------------------------------------------------------------------//
void IO_HDF5::write(SP_material mat){}

//-------------------------------------------------------------------------//
void IO_HDF5::write(SP_mesh mesh){}

//-------------------------------------------------------------------------//
void IO_HDF5::close(){}

//-------------------------------------------------------------------------//
IO_HDF5::SP_input IO_HDF5::read_input(){return SP_input(0);}

//-------------------------------------------------------------------------//
IO_HDF5::SP_material IO_HDF5::read_material(){return SP_material(0);}

//-------------------------------------------------------------------------//
IO_HDF5::SP_mesh IO_HDF5::read_mesh(){return SP_mesh(0);}

#endif // DETRAN_ENABLE_HDF5

} // end namespace detran_ioutils



//---------------------------------------------------------------------------//
//              end of file IO_HDF5.cc
//---------------------------------------------------------------------------//
