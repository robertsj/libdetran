#------------------------------------------------------------------------------#
# libdetran/src/solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran solvers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_SourceIteration test_SourceIteration.cc)
TARGET_LINK_LIBRARIES(test_SourceIteration 
                      solvers 
                      transport 
                      geometry 
                      materials 
                      utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

#ADD_TEST(test_SourceIteration_1D  test_SourceIteration 0)
ADD_TEST(test_SourceIteration_2D  test_SourceIteration 0)
#ADD_TEST(test_SourceIteration_3D  test_SourceIteration 3)