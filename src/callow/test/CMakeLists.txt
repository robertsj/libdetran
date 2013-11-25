#------------------------------------------------------------------------------#
# src/callow/test/CMakeLists.txt
#------------------------------------------------------------------------------#

# Vector
ADD_EXECUTABLE(test_Vector              test_Vector.cc)
TARGET_LINK_LIBRARIES(test_Vector       callow )
ADD_TEST(test_Vector                    test_Vector 0)

# Matrix
ADD_EXECUTABLE(test_Matrix              test_Matrix.cc)
TARGET_LINK_LIBRARIES(test_Matrix       callow )
ADD_TEST(test_Matrix                    test_Matrix 0)
#
ADD_EXECUTABLE(test_MatrixShell         test_MatrixShell.cc)
TARGET_LINK_LIBRARIES(test_MatrixShell  callow )
ADD_TEST(test_MatrixShell               test_MatrixShell 0)
#
ADD_EXECUTABLE(test_MatrixDense         test_MatrixDense.cc)
TARGET_LINK_LIBRARIES(test_MatrixDense  callow )
ADD_TEST(test_MatrixDense               test_MatrixDense 0)

# Linear Solvers
ADD_EXECUTABLE(test_LinearSolver        test_LinearSolver.cc)
TARGET_LINK_LIBRARIES(test_LinearSolver callow )
ADD_TEST(test_Richardson                test_LinearSolver 0)
ADD_TEST(test_Jacobi                    test_LinearSolver 1)
ADD_TEST(test_GaussSeidel               test_LinearSolver 2)
ADD_TEST(test_SOR                       test_LinearSolver 3)
ADD_TEST(test_GMRES                     test_LinearSolver 4)

# Eigenvalue Solvers
ADD_EXECUTABLE(test_EigenSolver         test_EigenSolver.cc)
TARGET_LINK_LIBRARIES(test_EigenSolver  callow )
ADD_TEST(test_PowerIteration            test_EigenSolver 0)

ADD_EXECUTABLE(test_Davidson            test_Davidson.cc)
TARGET_LINK_LIBRARIES(test_Davidson     callow )
ADD_TEST(test_Davidson_standard         test_Davidson   0)
ADD_TEST(test_Davidson_general          test_Davidson   0)

ADD_EXECUTABLE(test_Eispack             test_Eispack.cc)
TARGET_LINK_LIBRARIES(test_Eispack      callow )
ADD_TEST(test_Eispack                   test_Eispack 0)

# Preconditioners
ADD_EXECUTABLE(test_Preconditioners         test_Preconditioners.cc)
TARGET_LINK_LIBRARIES(test_Preconditioners  callow )
ADD_TEST(test_PCJacobi                      test_Preconditioners  0)
ADD_TEST(test_PCILU0                        test_Preconditioners  1)

# Performance, etc.
#ADD_EXECUTABLE(test_Threading               test_Threading.cc)
#TARGET_LINK_LIBRARIES(test_Threading        callow )
#ADD_TEST(test_ThreadedMatVec                test_Threading  0)
#ADD_TEST(test_ThreadedJacobi                test_Threading  1)

