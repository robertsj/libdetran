#------------------------------------------------------------------------------#
# libdetran/src/solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

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

# Test of Callow standard eigenvalue solvers (slepc potential)
ADD_EXECUTABLE(test_EigenArnoldi          test_EigenArnoldi.cc)
TARGET_LINK_LIBRARIES(test_EigenArnoldi   solvers)
ADD_TEST(test_EigenArnoldi_1g                  test_EigenArnoldi 0)
ADD_TEST(test_EigenArnoldi_7g_forward          test_EigenArnoldi 1)
ADD_TEST(test_EigenArnoldi_7g_adjoint          test_EigenArnoldi 2)

# Test of Callow's GD eigenvalue solver 
ADD_EXECUTABLE(test_EigenGD               test_EigenGD.cc)
TARGET_LINK_LIBRARIES(test_EigenGD   		  solvers)
ADD_TEST(test_EigenGD_1g                  test_EigenGD 0)
ADD_TEST(test_EigenGD_7g_forward          test_EigenGD 1)
ADD_TEST(test_EigenGD_7g_adjoint          test_EigenGD 2)

# Test of coarse mesh preconditioner
ADD_EXECUTABLE(test_MGCoarseMeshPreconditioner 
               test_MGCoarseMeshPreconditioner.cc)
TARGET_LINK_LIBRARIES(test_MGCoarseMeshPreconditioner solvers)
ADD_TEST(test_MGCoarseMeshPreconditioner_no_condense 
         test_MGCoarseMeshPreconditioner 0)
ADD_TEST(test_MGCoarseMeshPreconditioner_space
         test_MGCoarseMeshPreconditioner 1)
ADD_TEST(test_MGCoarseMeshPreconditioner_space_energy
         test_MGCoarseMeshPreconditioner 2)                  

# Test of MGDSA
ADD_EXECUTABLE(test_MGDSA test_MGDSA.cc)
TARGET_LINK_LIBRARIES(test_MGDSA solvers)
ADD_TEST(test_MGDSA_1g                  test_MGDSA 0)
ADD_TEST(test_MGDSA_7g_forward          test_MGDSA 1)
ADD_TEST(test_MGDSA_7g_forward_multiply test_MGDSA 2)
ADD_TEST(test_MGDSA_7g_adjoint          test_MGDSA 3)
ADD_TEST(test_MGDSA_7g_adjoint_multiply test_MGDSA 4)    

# Test of MGCMDSA
ADD_EXECUTABLE(test_MGCMDSA test_MGCMDSA.cc)
TARGET_LINK_LIBRARIES(test_MGCMDSA solvers)
ADD_TEST(test_MGCMDSA_1g                  test_MGCMDSA 0)
ADD_TEST(test_MGCMDSA_7g_forward          test_MGCMDSA 1)
ADD_TEST(test_MGCMDSA_7g_forward_multiply test_MGCMDSA 2)
ADD_TEST(test_MGCMDSA_7g_adjoint          test_MGCMDSA 3)
ADD_TEST(test_MGCMDSA_7g_adjoint_multiply test_MGCMDSA 4)    

# Test of SpectrumGS
ADD_EXECUTABLE(test_SpectrumGS test_SpectrumGS.cc)
TARGET_LINK_LIBRARIES(test_SpectrumGS solvers)
ADD_TEST(test_SpectrumGS                  test_SpectrumGS 0)
# Test of SpectrumPinCell
ADD_EXECUTABLE(test_SpectrumPinCell test_SpectrumPinCell.cc)
TARGET_LINK_LIBRARIES(test_SpectrumPinCell solvers)
ADD_TEST(test_SpectrumPinCell             test_SpectrumPinCell 0)


ADD_EXECUTABLE(test_MGSweepOperator          	test_MGSweepOperator.cc)
TARGET_LINK_LIBRARIES(test_MGSweepOperator   	solvers)
ADD_TEST(test_MGSweepOperator             		test_MGSweepOperator 0)

ADD_EXECUTABLE(test_MGCMFD          	test_MGCMFD.cc)
TARGET_LINK_LIBRARIES(test_MGCMFD   	solvers)
ADD_TEST(test_MGCMFD             		test_MGCMFD 0)


ADD_EXECUTABLE(test_TimeStepper          	test_TimeStepper.cc)
TARGET_LINK_LIBRARIES(test_TimeStepper   	solvers)
ADD_EXECUTABLE(test_TWIGL                   test_TWIGL.cc)
TARGET_LINK_LIBRARIES(test_TWIGL            solvers)
#ADD_TEST(test_MGCMFD             		test_MGCMFD 0)
