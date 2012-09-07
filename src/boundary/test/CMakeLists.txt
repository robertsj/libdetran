#------------------------------------------------------------------------------#
# src/boundary/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran boundary tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_BoundarySN              test_BoundarySN.cc)
TARGET_LINK_LIBRARIES(BoundarySN            boundary)

ADD_EXECUTABLE(test_BoundaryMOC             test_BoundaryMOC.cc)
TARGET_LINK_LIBRARIES(BoundaryMOC           boundary)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_BoundarySN    test_BoundarySN       0)
ADD_TEST(test_BoundaryMOC   test_BoundaryMOC      0)

