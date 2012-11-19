#------------------------------------------------------------------------------#
# src/discretization/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran discretization tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Equation_DD_1D              test_Equation_DD_1D.cc)
TARGET_LINK_LIBRARIES(test_Equation_DD_1D       discretization)
ADD_EXECUTABLE(test_Equation_SC_1D              test_Equation_SC_1D.cc)
TARGET_LINK_LIBRARIES(test_Equation_SC_1D       discretization)

#ADD_EXECUTABLE(test_BoundaryMOC                 test_BoundaryMOC.cc)
#TARGET_LINK_LIBRARIES(test_BoundaryMOC          boundary)

#ADD_EXECUTABLE(test_BoundaryDiffusion           test_BoundaryDiffusion.cc)
#TARGET_LINK_LIBRARIES(test_BoundaryDiffusion    boundary)

##------------------------------------------------------------------------------#
## Define tests
##------------------------------------------------------------------------------#

ADD_TEST(test_Equation_DD_1D          test_Equation_DD_1D         0)
ADD_TEST(test_Equation_SC_1D          test_Equation_SC_1D         0)
#ADD_TEST(test_BoundaryMOC           test_BoundaryMOC        0)
#ADD_TEST(test_BoundaryDiffusion     test_BoundaryDiffusion  0)
