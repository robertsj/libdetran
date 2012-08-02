#------------------------------------------------------------------------------#
# libdetran/src/ioutils/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran ioutils tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

if (DETRAN_ENABLE_HDF5)

ADD_EXECUTABLE(test_IO_HDF5 test_IO_HDF5.cc)
TARGET_LINK_LIBRARIES(test_IO_HDF5 ioutils)

endif (DETRAN_ENABLE_HDF5)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_IO_HDF5_input       test_IO_HDF5   0)
ADD_TEST( test_IO_HDF5_material    test_IO_HDF5   1)
