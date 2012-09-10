#------------------------------------------------------------------------------#
# libdetran/src/solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran solvers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_SourceIteration 			test_SourceIteration.cc)
TARGET_LINK_LIBRARIES(test_SourceIteration 		solvers )

ADD_EXECUTABLE(test_PowerIteration 				test_PowerIteration.cc)
TARGET_LINK_LIBRARIES(test_PowerIteration  		solvers)

if (DETRAN_ENABLE_PETSC)

ADD_EXECUTABLE(test_InnerGMRES 					test_InnerGMRES.cc)
TARGET_LINK_LIBRARIES(test_InnerGMRES 			solvers)

ADD_EXECUTABLE(test_OperatorMatrix 			    test_OperatorMatrix.cc)
TARGET_LINK_LIBRARIES(test_OperatorMatrix		solvers)

ADD_EXECUTABLE(test_OperatorShell			    test_OperatorShell.cc)
TARGET_LINK_LIBRARIES(test_OperatorShell		solvers)

endif ()

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

#ADD_TEST(test_SourceIteration_1D  test_SourceIteration 0)
ADD_TEST(test_SourceIteration_2D  test_SourceIteration 0)
#ADD_TEST(test_SourceIteration_3D  test_SourceIteration 3)
ADD_TEST(test_PowerIteration_2D test_PowerIteration 0)

if (DETRAN_ENABLE_PETSC)

ADD_TEST(test_InnerGMRES_1D 		test_InnerGMRES 0)
ADD_TEST(test_OperatorMatrix 		test_OperatorMatrix 0)
ADD_TEST(test_OperatorShell 		test_OperatorShell 0)

endif ()
