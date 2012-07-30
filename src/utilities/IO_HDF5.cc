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
{

}

void IO_HDF5::write(SP_input input)
{
  // Preconditions
  Require(input);

  // Open the HDF5 file
  d_file_id = H5Fcreate(d_filename,     // filename
                        H5F_ACC_TRUNC,  // overwrite existing file
                        H5P_DEFAULT,    // file create property list
                        H5P_DEFAULT);   // file access property list




  // Close the HDF5 file
  herr_t status = H5Fclose(d_file_id);
  if (!status) std::cout << "Error closing HDF5 file"+d_filename << std::endl;
}

void IO_HDF5::close()
{

}

void IO_HDF5::read(SP_input input, std::string filename)
{

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file IO_HDF5.cc
//---------------------------------------------------------------------------//
