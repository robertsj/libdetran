//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   IO_HDF5.t.hh
 *  @author robertsj
 *  @date   Aug 1, 2012
 *  @brief  IO_HDF5 template definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_ioutils_IO_HDF5_T_HH_
#define detran_ioutils_IO_HDF5_T_HH_

namespace detran_ioutils
{

//---------------------------------------------------------------------------//
template <class T>
hid_t IO_HDF5::set_memtype()
{
  // Define the compound data type.
  typedef compound_type<T> data_type;

  // Define the memory type for the compound type's value.
  HDF5_MemoryType mem;
  hid_t value_memtype = mem.type<T>();

  // Define the memory type for the compound type's key.
  HDF5_MemoryType stringmem;
  hid_t string_memtype = stringmem.type<std::string>();

  // Define the compound type's memory type
  hid_t memtype;
  herr_t status;
  memtype = H5Tcreate(H5T_COMPOUND, sizeof(data_type));
  status  = H5Tinsert(memtype, "KEY",
                      HOFFSET(data_type, key),   string_memtype);
  status  = H5Tinsert(memtype, "VALUE",
                      HOFFSET(data_type, value), value_memtype);
  return memtype;
}

//---------------------------------------------------------------------------//
template <class T>
hid_t IO_HDF5::set_filetype()
{
  // Define the compound data type.
  typedef compound_type<T> data_type;

  // Define the memory type for the compound type's value.
  HDF5_FileType mem;
  hid_t value_memtype = mem.type<T>();

  // Define the memory type for the compound type's key.
  HDF5_FileType stringmem;
  hid_t string_memtype = stringmem.type<std::string>();

  // Compound type size
  size_t size = sizeof (hvl_t) +
                H5Tget_size(value_memtype);

  // Define the compound type's memory type
  hid_t filetype;
  herr_t status;
  filetype = H5Tcreate (H5T_COMPOUND, size);
  status = H5Tinsert (filetype, "KEY",   0,              string_memtype);
  status = H5Tinsert (filetype, "VALUE", sizeof (hvl_t), value_memtype);
  return filetype;
}

//---------------------------------------------------------------------------//
template <class T>
bool IO_HDF5::read_data(SP_input 		input,
                        hid_t 			group,
                        std::string name)
{
  // Preconditions
  Require(input);

  // Open the dataset
  hid_t dset = H5Dopen(group, name.c_str(), H5P_DEFAULT);

  // Get dataspace and allocate memory for read buffer.
  hsize_t dims[1];
  hid_t space = H5Dget_space (dset);
  int ndims   = H5Sget_simple_extent_dims(space, dims, NULL);
  herr_t status;

  // Return if there is no data.
  if (!dims[0])
  {
    status = H5Sclose(space);
    return false;
  }

  // Create buffer for reading.
  compound_type<T>* data(new compound_type<T>[dims[0]]);

  // Create the compound datatype for memory.
  hid_t memtype = set_memtype<T>();

  // Read the data.
  status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Loop and insert into input
  for (int i = 0; i < dims[0]; i++)
  {
    std::string key(data[i].key);
    T value;
    get_value(value, data[i].value);
    input->put(key, value);
  }

  // Close and release resources.  H5Dvlen_reclaim will automatically
  // traverse the structure and free any vlen data.
  status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, data);
  delete [] data;
  status = H5Dclose(dset);
  status = H5Sclose(space);

  return true;
}

//---------------------------------------------------------------------------//
template <>
bool IO_HDF5::read_data<IO_HDF5::SP_input>(SP_input      db,
                                           hid_t         root,
                                           std::string   name);

//---------------------------------------------------------------------------//
template <class T>
bool IO_HDF5::write_map(hid_t group,
                        const char *name,
                        const std::map<std::string, T> &map)
{
  // Preconditions
  /* ... */

  // Get the size of this data set
  int size = map.size();

  // Create a buffer for writing out.
  compound_type<T>* data(new compound_type<T>[size]);

  // Loop and add
  typename std::map<std::string, T>::const_iterator it = map.begin();
  int i = 0;
  for (; it != map.end(); it++, i++)
  {
    set_value(it->first,  data[i].key);
    set_value(it->second, data[i].value);
  }
  Assert(i == size);

  // Set dimension
  hsize_t dims[1] = {size};

  // Create the compound datatype for memory.
  hid_t memtype = set_memtype<T>();

  // Create the compound datatype for file.
  hid_t filetype = set_filetype<T>();

  // Create the dataspace, the dataset, and write to file.
  hid_t space   = H5Screate_simple(1, dims, NULL);
  hid_t dset    = H5Dcreate(group, name, filetype, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Only write data if there is some.  This lets an empty data set remain.
  herr_t status;
  if (size > 0)
  {
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }

  // Close out everything
  delete [] data;
  status = H5Dclose(dset);
  status = H5Sclose(space);
  status = H5Tclose(filetype);
  status = H5Tclose(memtype);

  return true;
}

//---------------------------------------------------------------------------//
template <>
inline bool IO_HDF5::write_map<IO_HDF5::SP_input>
   (hid_t rootgroup,
		const char *name,
		const std::map<std::string, SP_input> &map)
{
  // Create the db group
  hid_t group = H5Gcreate(rootgroup, name,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Get the size of this data set
  int size = map.size();

  // Loop and add
  typedef std::map<std::string, SP_input> map_T;
  map_T::const_iterator it = map.begin();
  int i = 0;
  for (; it != map.end(); ++it++, ++i)
  {
  	write(it->second, it->first, group);
  }
  Assert(i == size);

  return true;
}

//---------------------------------------------------------------------------//
template <class T>
bool IO_HDF5::read_map(hid_t 										 group,
                       const char 							*name,
                       std::map<std::string, T> &map)
{

  // Ensure dataset is present.
  if (!exists(group, name)) return false;

  // Open the dataset
  hid_t dset = H5Dopen(group, name, H5P_DEFAULT);

  // Get dataspace and allocate memory for read buffer.
  hsize_t dims[1];
  hid_t space = H5Dget_space (dset);
  int ndims   = H5Sget_simple_extent_dims(space, dims, NULL);

  // Return if there is no data.
  if (!dims[0]) return false;

  // Create buffer for reading.
  compound_type<T>* data(new compound_type<T>[dims[0]]);

  // Create the compound datatype for memory.
  hid_t memtype = set_memtype<T>();

  // Read the data.
  herr_t status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Loop and insert into input
  for (int i = 0; i < dims[0]; i++)
  {
    std::string key(data[i].key);
    T value;
    get_value(value, data[i].value);
    map[key] = value;
  }

  // Close and release resources.  H5Dvlen_reclaim will automatically
  // traverse the structure and free any vlen data.
  status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, data);
  delete [] data;
  status = H5Dclose(dset);
  status = H5Sclose(space);


  return true;
}

//---------------------------------------------------------------------------//
template <class T>
bool IO_HDF5::write_vec(hid_t group, const char* name, const std::vector<T> &source)
{
  // Preconditions
  Require(source.size());

  // Datasource dimension.  This is only for 1D data.
  hsize_t dims[1] = {source.size()};

  // Create the dataspace.
  hid_t space = H5Screate_simple(1, dims, NULL);

  // Create the dataset, write to it, and close it.
  HDF5_MemoryType mem;
  hid_t dset = H5Dcreate(group, name, mem.type<T>(), space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, &source[0]);
  status = H5Dclose(dset);

  return true;
}

//---------------------------------------------------------------------------//
template <class T>
bool IO_HDF5::read_vec(hid_t group, const char* name, std::vector<T> &target)
{
  // Preconditions
  Require(target.size());

  // Create buffer for reading.
  T* buffer(new T[target.size()]);

  // Ensure dataset is present.
  if (!exists(group, name)) return false;

  HDF5_MemoryType mem;
  hid_t dset = H5Dopen(group, name, H5P_DEFAULT);
  herr_t status = H5Dread(dset, mem.type<T>(), H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, buffer);
  status = H5Dclose(dset);

  // Copy buffer to target vector.
  for (int i = 0; i < target.size(); i++)
    target[i] = buffer[i];

  // Delete buffer
  delete [] buffer;

  return true;
}

//---------------------------------------------------------------------------//
template <class T>
bool IO_HDF5::write_scalar_attribute
(hid_t group, const char* name, const T &value)
{
  // Create scalar attribute.
  HDF5_MemoryType mem;
  hid_t att_space = H5Screate(H5S_SCALAR);
  hid_t att = H5Acreate(group, name, mem.type<T>(), att_space,
                        H5P_DEFAULT, H5P_DEFAULT);

  // Write scalar attribute and close everything.
  herr_t status = H5Awrite(att, H5T_NATIVE_INT, &value);
  status = H5Sclose(att_space);
  status = H5Aclose(att);

  return true;
}

//---------------------------------------------------------------------------//
template <class T>
bool  IO_HDF5::read_scalar_attribute
(hid_t group, const char* name, T &value)
{
  hid_t att = H5Aopen_name(group, name);
  herr_t status = H5Aread(att, H5T_NATIVE_INT, &value);
  status = H5Aclose(att);
  return true;
}

} // end namespace detran_ioutils

#endif /* IO_HDF5_T_HH_ */
