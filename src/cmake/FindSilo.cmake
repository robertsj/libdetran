# (Slightly adapted from S. Johnson)
# - Find LLNL's Silo library
# This module defines
#   Silo_INCLUDE_DIR, where to find blitz/blitz.h, etc.
#   Silo_LIBRARIES, libraries to link against to use Silo.
#   Silo_FOUND, If false, do not try to use Silo.
# also defined, but not for general use are
#   Silo_LIBRARY, where to find the Silo library.
# The user should specify the head Silo director, Silo_DIR,
# or Silo_INC and Silo_LIB.  

if (Silo_DIR)
    message( "-- Using Silo_DIR: ${Silo_DIR} " )
endif()
if (Silo_LIB)
    message( "-- Using Silo_LIB: ${Silo_LIB} " )
endif()
if (Silo_INC)
    message( "-- Using Silo_INC: ${Silo_INC} " )
endif()

# Set the cmake library path.  This is because a clean
# build of cmake 2.8.4 seems to exclude this path.  It might
# be a Ubuntu-specific quirk.
set(CMAKE_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")

find_path(Silo_INCLUDE_DIR 
          NAMES silo.h
          PATHS ${Silo_DIR}/include 
                ${Silo_INC}
)
message("-- Using Silo_INCLUDE_DIR: ${Silo_INCLUDE_DIR} " )

find_library(Silo_LIBRARY
             NAMES siloh5 silo siloxx
             PATHS ${Silo_DIR}/lib
                   ${Silo_LIB}
)
message("-- Using Silo_LIBRARY: ${Silo_LIBRARY} " )

if (Silo_LIBRARY MATCHES "siloh5")
    MESSAGE( "-- Note: This Silo installation needs HDF5." )
    FIND_PACKAGE(HDF5 REQUIRED)
endif()

# Requires ZLib
find_package(ZLIB REQUIRED)

SET( Silo_FOUND "NO" )
IF(Silo_INCLUDE_DIR)
  IF(Silo_LIBRARY)

    SET( Silo_LIBRARIES ${Silo_LIBRARY})
    SET( Silo_FOUND "YES" )

    #The following deprecated settings are for backwards compatibility with CMake1.4
    SET (Silo_INCLUDE_PATH ${Silo_INCLUDE_DIR})

  ENDIF(Silo_LIBRARY)
ENDIF(Silo_INCLUDE_DIR)

MARK_AS_ADVANCED(
  Silo_INCLUDE_DIR
  Silo_LIBRARY
)


