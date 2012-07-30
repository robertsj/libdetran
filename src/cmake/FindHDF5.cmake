#
# Find the native HDF5 includes and library
#
# HDF5_INCLUDE_DIR - where to find H5public.h, etc.
# HDF5_LIBRARIES   - List of fully qualified libraries to link against when using hdf5.
# HDF5_FOUND       - Do not attempt to use hdf5 if "no" or undefined.

find_path(HDF5_INCLUDE_DIR 
          NAMES H5public.h
          PATHS ${HDF5_DIR}/include 
                ${HDF5_INC}
)
message("HDF5 Include: ${HDF5_INCLUDE_DIR}")

find_library(HDF5_LIBRARY hdf5
             PATHS ${HDF5_DIR}/lib
                   ${HDF5_LIB}
)
message("HDF5 Library: ${HDF5_LIBRARY}")

set( HDF5_FOUND "NO" )
if(HDF5_INCLUDE_DIR)
  if(HDF5_LIBRARY)
    set( HDF5_LIBRARIES ${HDF5_LIBRARY})
    set( HDF5_FOUND "YES" )
  endif(HDF5_LIBRARY)
elseif(HDF5_INCLUDE_DIR)
  message("HDF5 not found. Try setting HDF5_DIR.")
endif(HDF5_INCLUDE_DIR)


