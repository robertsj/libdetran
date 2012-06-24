#------------------------------------------------------------------------------#
# libdetran/src/geometry/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran geometry tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Mesh1D test_Mesh1D.cc)
TARGET_LINK_LIBRARIES(test_Mesh1D geometry utilities)

ADD_EXECUTABLE(test_Mesh2D test_Mesh2D.cc)
TARGET_LINK_LIBRARIES(test_Mesh2D geometry utilities)

ADD_EXECUTABLE(test_Mesh3D test_Mesh3D.cc)
TARGET_LINK_LIBRARIES(test_Mesh3D geometry utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_Mesh1D 	test_Mesh1D 	0)
ADD_TEST(test_Mesh2D 	test_Mesh2D 	0)
ADD_TEST(test_Mesh3D 	test_Mesh3D 	0)
