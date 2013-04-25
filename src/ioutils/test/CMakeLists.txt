#------------------------------------------------------------------------------#
# libdetran/src/ioutils/test/CMakeLists.txt
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

if (DETRAN_ENABLE_HDF5)
ADD_EXECUTABLE(test_IO_HDF5 test_IO_HDF5.cc)
TARGET_LINK_LIBRARIES(test_IO_HDF5 ioutils)
endif ()

if(DETRAN_ENABLE_SILO)
ADD_EXECUTABLE(test_SiloOutput test_SiloOutput.cc)
TARGET_LINK_LIBRARIES(test_SiloOutput ioutils)
endif ()

ADD_EXECUTABLE(test_PPMPlotter test_PPMPlotter.cc)
TARGET_LINK_LIBRARIES(test_PPMPlotter ioutils)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

if (DETRAN_ENABLE_HDF5)
ADD_TEST( test_IO_HDF5_input       test_IO_HDF5   0)
ADD_TEST( test_IO_HDF5_material    test_IO_HDF5   1)
ADD_TEST( test_IO_HDF5_mesh        test_IO_HDF5   2)
endif ()

if(DETRAN_ENABLE_SILO)
ADD_TEST( test_SiloOutput  test_SiloOutput 0)
endif()
