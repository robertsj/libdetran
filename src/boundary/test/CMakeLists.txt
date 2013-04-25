#------------------------------------------------------------------------------#
# src/boundary/test/CMakeLists.txt
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_BoundarySN                  test_BoundarySN.cc)
TARGET_LINK_LIBRARIES(test_BoundarySN           boundary)

ADD_EXECUTABLE(test_BoundaryMOC                 test_BoundaryMOC.cc)
TARGET_LINK_LIBRARIES(test_BoundaryMOC          boundary)

ADD_EXECUTABLE(test_BoundaryDiffusion           test_BoundaryDiffusion.cc)
TARGET_LINK_LIBRARIES(test_BoundaryDiffusion    boundary)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_BoundarySN            test_BoundarySN         0)
ADD_TEST(test_BoundaryMOC           test_BoundaryMOC        0)
ADD_TEST(test_BoundaryDiffusion     test_BoundaryDiffusion  0)
