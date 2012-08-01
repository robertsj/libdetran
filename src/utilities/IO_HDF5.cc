//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   IO_HDF5.cc
 * \brief  IO_HDF5 
 * \author Jeremy Roberts
 * \date   Jul 29, 2012
 */
//---------------------------------------------------------------------------//

#include "IO_HDF5.hh"

namespace detran
{

IO_HDF5::IO_HDF5(std::string filename)
  : d_filename(filename)
  , d_open(false)
{
  /* ... */
}

void IO_HDF5::write(SP_input input)
{
  // Preconditions
  Require(input);

  // Open the HDF5 file
  d_file_id = H5Fcreate(d_filename.c_str(), // filename
                        H5F_ACC_TRUNC,      // overwrite existing file
                        H5P_DEFAULT,        // file create property list
                        H5P_DEFAULT);       // file access property list
  d_open = true;

  // Create the input group
  hid_t group = H5Gcreate(d_file_id, "/input",
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t  status;
  hid_t memtype;
  hid_t filetype;
  hid_t dset;
  hid_t space;

  // Create variable string type
  d_string_type = H5Tcopy (H5T_C_S1);
  status = H5Tset_size (d_string_type, H5T_VARIABLE);

  //-------------------------------------------------------------------------//
  // PREPARE DATA
  //-------------------------------------------------------------------------//

  // Get the map sizes for the suppored types.
  int size_int      = input->size(InputDB::INT);
  int size_dbl      = input->size(InputDB::DBL);
  int size_str      = input->size(InputDB::STR);
  int size_vec_int  = input->size(InputDB::VEC_INT);
  int size_vec_dbl  = input->size(InputDB::VEC_DBL);

  // Define the temporary data containers.
  compound_type<int>            data_int[size_int];
  compound_type<double>         data_dbl[size_dbl];
  compound_type<std::string>    data_str[size_str];
  compound_type<vec_int>        data_vec_int[size_vec_int];
  compound_type<vec_dbl>        data_vec_dbl[size_vec_dbl];

  //-------------------------------------------------------------------------//
  // WRITE DATA
  //-------------------------------------------------------------------------//

  write_data(input, group, "int_data", data_int);
  write_data(input, group, "dbl_data", data_dbl);
  write_data(input, group, "str_data", data_str);
  write_data(input, group, "vec_int_data", data_vec_int);
  write_data(input, group, "vec_dbl_data", data_vec_dbl);

  // Close the group.
  status = H5Gclose(group);

}

void IO_HDF5::close()
{
  // If not open, ignore.
  if (!d_open) return;

  // Close the HDF5 file
  H5Fclose(d_file_id);
  d_open = false;
}

void IO_HDF5::read(SP_input input)
{

  // Open the file if necessary.
  if (!d_open)
    d_file_id = H5Fopen(d_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the input group
  hid_t group = H5Gopen(d_file_id, "input", H5P_DEFAULT);

  // Read data
  read_data<int>(input, group, "int_data");
  read_data<double>(input, group, "dbl_data");
  read_data<std::string>(input, group, "str_data");
  read_data<vec_int>(input, group, "vec_int_data");
  read_data<vec_dbl>(input, group, "vec_dbl_data");

  // Close the group.
  herr_t status = H5Gclose(group);

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file IO_HDF5.cc
//---------------------------------------------------------------------------//
