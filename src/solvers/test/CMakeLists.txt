#------------------------------------------------------------------------------#
# libdetran/src/solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_FixedSourceManager          test_FixedSourceManager.cc)
TARGET_LINK_LIBRARIES(test_FixedSourceManager   solvers)

ADD_EXECUTABLE(test_EigenvalueManager           test_EigenvalueManager.cc)
TARGET_LINK_LIBRARIES(test_EigenvalueManager    solvers)

ADD_EXECUTABLE(test_TimeStepper           		test_TimeStepper.cc)
TARGET_LINK_LIBRARIES(test_TimeStepper    		solvers)

ADD_EXECUTABLE(test_TWIGL                      	test_TWIGL.cc)
TARGET_LINK_LIBRARIES(test_TWIGL              	solvers)

ADD_EXECUTABLE(test_LRA                      	test_LRA.cc)
TARGET_LINK_LIBRARIES(test_LRA              	solvers)

# Test of Gauss-Seidel
ADD_EXECUTABLE(test_MGSolverGS                  test_MGSolverGS.cc)
TARGET_LINK_LIBRARIES(test_MGSolverGS           solvers)
ADD_TEST(test_MGSolverGS_1g                     test_MGSolverGS 0)
ADD_TEST(test_MGSolverGS_7g_forward             test_MGSolverGS 1)
ADD_TEST(test_MGSolverGS_7g_forward_multiply    test_MGSolverGS 2)
ADD_TEST(test_MGSolverGS_7g_adjoint             test_MGSolverGS 3)
ADD_TEST(test_MGSolverGS_7g_adjoint_multiply    test_MGSolverGS 4)

# Test of Multigroup GMRES
ADD_EXECUTABLE(test_MGSolverGMRES               test_MGSolverGMRES.cc)
TARGET_LINK_LIBRARIES(test_MGSolverGMRES        solvers)
ADD_TEST(test_MGSolverGMRES_1g                  test_MGSolverGMRES 0)
ADD_TEST(test_MGSolverGMRES_7g_forward          test_MGSolverGMRES 1)
ADD_TEST(test_MGSolverGMRES_7g_forward_multiply test_MGSolverGMRES 2)
ADD_TEST(test_MGSolverGMRES_7g_adjoint          test_MGSolverGMRES 3)
ADD_TEST(test_MGSolverGMRES_7g_adjoint_multiply test_MGSolverGMRES 4)

# Test of Multigroup Diffusion
ADD_EXECUTABLE(test_MGDiffusionSolver               test_MGDiffusionSolver.cc)
TARGET_LINK_LIBRARIES(test_MGDiffusionSolver        solvers)
ADD_TEST(test_MGDiffusionSolver_1g                  test_MGDiffusionSolver 0)
ADD_TEST(test_MGDiffusionSolver_7g_forward          test_MGDiffusionSolver 1)
ADD_TEST(test_MGDiffusionSolver_7g_forward_multiply test_MGDiffusionSolver 2)
ADD_TEST(test_MGDiffusionSolver_7g_adjoint          test_MGDiffusionSolver 3)
ADD_TEST(test_MGDiffusionSolver_7g_adjoint_multiply test_MGDiffusionSolver 4)

# Test of Power Iteration
ADD_EXECUTABLE(test_EigenPI               test_EigenPI.cc)
TARGET_LINK_LIBRARIES(test_EigenPI        solvers)
ADD_TEST(test_EigenPI_1g                  test_EigenPI 0)
ADD_TEST(test_EigenPI_7g_forward          test_EigenPI 1)
ADD_TEST(test_EigenPI_7g_adjoint          test_EigenPI 2)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

#ADD_TEST(test_SourceIteration_1D  test_SourceIteration 0)
#ADD_TEST(test_SourceIteration_2D  test_SourceIteration 0)
#ADD_TEST(test_SourceIteration_3D  test_SourceIteration 3)
#ADD_TEST(test_PowerIteration_2D test_PowerIteration 0)
ADD_TEST(test_FixedSourceManager_1D        test_FixedSourceManager 0)
ADD_TEST(test_EigenvalueManager_1D         test_EigenvalueManager 0)
ADD_TEST(test_TimeStepper_BDF              test_TimeStepper 1)
#ADD_TEST(test_TWIGL              test_TWIGL 0)
