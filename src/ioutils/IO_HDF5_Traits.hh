//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   IO_HDF5_Traits.hh
 * \author robertsj
 * \date   Aug 1, 2012
 * \brief  IO_HDF5_Traits class definition.
 * \todo   It may be good to extract the classes and such.
 */
//---------------------------------------------------------------------------//

#ifndef IO_HDF5_TRAITS_HH_
#define IO_HDF5_TRAITS_HH_

// Utilities
#include "DBC.hh"

// System
#include <iostream>
#include <string>
#include <vector>
#include "hdf5.h"

namespace detran_ioutils
{

//---------------------------------------------------------------------------//

/// Compound type traits
template <class T>
struct compound_type_traits
{
  typedef T value_type;
};
template <>
struct compound_type_traits<std::string>
{
  typedef char* value_type;
};
template <>
struct compound_type_traits<std::vector<int> >
{
  typedef hvl_t value_type;
};
template <>
struct compound_type_traits<std::vector<double> >
{
  typedef hvl_t value_type;
};

//---------------------------------------------------------------------------//

/// Compound type
template <class T>
struct compound_type
{
  char* key;
  typename compound_type_traits<T>::value_type value;
};

//---------------------------------------------------------------------------//

/*!
 *  \brief Set the value
 *
 *  Given the actual data type as stored in the map, give the new
 *  data for storage.
 */
template <class T>
void set_value(const T &value_in,
               typename compound_type_traits<T>::value_type &value_out)
{
  value_out = value_in;
}
template <>
void set_value(const std::string &value_in,
               compound_type_traits<std::string>::value_type &value_out)
{
  value_out = (char*) value_in.c_str();
}
template <>
void set_value(const std::vector<int> &value_in,
               compound_type_traits<std::vector<int> >::value_type &value_out)
{
  value_out.len = value_in.size();
  value_out.p   = (void *) &value_in[0];
}
template <>
void set_value(const std::vector<double> &value_in,
               compound_type_traits<std::vector<double> >::value_type &value_out)
{
  value_out.len = value_in.size();
  value_out.p   = (void *) &value_in[0];
}

/*!
 *  \brief Get the value
 */
template <class T>
void get_value(T &value_out,
               const typename compound_type_traits<T>::value_type &value_in)
{
  value_out = value_in;
}
template <>
void get_value(std::string &value_out,
               const compound_type_traits<std::string>::value_type &value_in)
{
  value_out = (char*) value_in;
}
template <>
void get_value(std::vector<int> &value_out,
               const compound_type_traits<std::vector<int> >::value_type &value_in)
{
  // Resize the vector.
  value_out.resize(value_in.len);

  // Get an integer pointer to the input data.
  int *tmp = (int *) value_in.p;

  // Fill the vector
  for (int i = 0; i < value_out.size(); i++)
  {
    value_out[i] = tmp[i];
  }
}
template <>
void get_value(std::vector<double> &value_out,
               const compound_type_traits<std::vector<double> >::value_type &value_in)
{
  // Resize the vector.
  value_out.resize(value_in.len);

  // Get an integer pointer to the input data.
  double *tmp = (double *) value_in.p;

  // Fill the vector
  for (int i = 0; i < value_out.size(); i++)
  {
    value_out[i] = tmp[i];
  }
}


//---------------------------------------------------------------------------//
/*!
 *  \class HDF5_MemoryType
 *  \brief Create the memory type for the map entry value as stored in memory.
 */
//---------------------------------------------------------------------------//
class HDF5_MemoryType
{

public:

  /// Destructor.  Note, the type must be closed to avoid a leak.
  ~HDF5_MemoryType()
  {
    H5Tclose(d_type);
  }

  /// Return the type for the map entry value.
  template <class T>
  hid_t type();

private:

  /// HDF5 type for the map entry value.
  hid_t d_type;
};

template <>
hid_t HDF5_MemoryType::type<int>()
{
  d_type = H5Tcopy(H5T_NATIVE_INT);
  return d_type;
}
template <>
hid_t HDF5_MemoryType::type<double>()
{
  d_type = H5Tcopy(H5T_NATIVE_DOUBLE);
  return d_type;
}
template <>
hid_t HDF5_MemoryType::type<std::string>()
{
  d_type = H5Tcopy(H5T_C_S1);
  herr_t status = H5Tset_size (d_type, H5T_VARIABLE);
  return d_type;
}
template <>
hid_t HDF5_MemoryType::type<std::vector<int> >()
{
  d_type = H5Tvlen_create(H5T_NATIVE_INT);
  return d_type;
}
template <>
hid_t HDF5_MemoryType::type<std::vector<double> >()
{
  d_type = H5Tvlen_create(H5T_NATIVE_DOUBLE);
  return d_type;
}

//---------------------------------------------------------------------------//
/*!
 *  \class HDF5_MemoryType
 *  \brief Create the memory type for the map entry value as stored on disk.
 */
//---------------------------------------------------------------------------//
class HDF5_FileType
{

public:

  /// Destructor.  Note, the type must be closed to avoid a leak.
  ~HDF5_FileType(){ H5Tclose(d_type); }

  /// Return the type for the map entry value.
  template <class T>
  hid_t type();

private:

  /// HDF5 type for the map entry value.
  hid_t d_type;
};

template <>
hid_t HDF5_FileType::type<int>()
{
  d_type = H5Tcopy(H5T_STD_I32LE);
  return d_type;
}
template <>
hid_t HDF5_FileType::type<double>()
{
  d_type = H5Tcopy(H5T_IEEE_F64BE);
  return d_type;
}
template <>
hid_t HDF5_FileType::type<std::string>()
{
  d_type = H5Tcopy(H5T_C_S1);
  herr_t status = H5Tset_size (d_type, H5T_VARIABLE);
  return d_type;
}
template <>
hid_t HDF5_FileType::type<std::vector<int> >()
{
  d_type = H5Tvlen_create(H5T_STD_I32LE);
  return d_type;
}
template <>
hid_t HDF5_FileType::type<std::vector<double> >()
{
  d_type = H5Tvlen_create(H5T_IEEE_F64BE);
  return d_type;
}

//



} // end namespace detran_ioutils

#endif /* IO_HDF5_TRAITS_HH_ */
