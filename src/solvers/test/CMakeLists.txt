#------------------------------------------------------------------------------#
# libdetran/src/solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran solvers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

#ADD_EXECUTABLE(test_SourceIteration 			test_SourceIteration.cc)
#TARGET_LINK_LIBRARIES(test_SourceIteration 		solvers )

#ADD_EXECUTABLE(test_PowerIteration 				test_PowerIteration.cc)
#TARGET_LINK_LIBRARIES(test_PowerIteration  		solvers)


#ADD_EXECUTABLE(test_InnerGMRES                 test_InnerGMRES.cc)
#TARGET_LINK_LIBRARIES(test_InnerGMRES          solvers)

#ADD_EXECUTABLE(test_DiffusionFixedSourceSolver          test_DiffusionFixedSourceSolver.cc)
#TARGET_LINK_LIBRARIES(test_DiffusionFixedSourceSolver   solvers)

#ADD_EXECUTABLE(test_DiffusionEigensolver          test_DiffusionEigensolver.cc)
#TARGET_LINK_LIBRARIES(test_DiffusionEigensolver   solvers)

ADD_EXECUTABLE(test_FixedSourceManager          test_FixedSourceManager.cc)
TARGET_LINK_LIBRARIES(test_FixedSourceManager   solvers)

ADD_EXECUTABLE(test_EigenvalueManager           test_EigenvalueManager.cc)
TARGET_LINK_LIBRARIES(test_EigenvalueManager    solvers)


#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

#ADD_TEST(test_SourceIteration_1D  test_SourceIteration 0)
#ADD_TEST(test_SourceIteration_2D  test_SourceIteration 0)
#ADD_TEST(test_SourceIteration_3D  test_SourceIteration 3)
#ADD_TEST(test_PowerIteration_2D test_PowerIteration 0)
ADD_TEST(test_FixedSourceManager_1D        test_FixedSourceManager 0)
ADD_TEST(test_EigenvalueManager_1D         test_EigenvalueManager 0)

if (DETRAN_ENABLE_PETSC)

#ADD_TEST(test_InnerGMRES_1D 		test_InnerGMRES 0)
#ADD_TEST(test_DiffusionFixedSourceSolver_1D  test_DiffusionFixedSourceSolver 0)
#ADD_TEST(test_DiffusionFixedSourceSolver_2D  test_DiffusionFixedSourceSolver 1)
#ADD_TEST(test_DiffusionFixedSourceSolver_3D  test_DiffusionFixedSourceSolver 2)
#ADD_TEST(test_DiffusionEigenSolver_1D        test_DiffusionEigensolver 0)

endif ()
