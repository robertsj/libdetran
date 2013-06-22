#------------------------------------------------------------------------------#
# libdetran/src/utilities/test/CMakeLists.txt
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Testing                 test_Testing.cc)
TARGET_LINK_LIBRARIES(test_Testing          utilities)
ADD_TEST(test_Testing_pass                  test_Testing 0)
ADD_TEST(test_Testing_fail                  test_Testing 1)


ADD_EXECUTABLE(test_InputDB                 test_InputDB.cc)
TARGET_LINK_LIBRARIES(test_InputDB          utilities)
ADD_TEST(test_InputDB_basic                 test_InputDB 0)

ADD_EXECUTABLE(test_MathUtilities           test_MathUtilities.cc)
TARGET_LINK_LIBRARIES(test_MathUtilities    utilities)
ADD_TEST( test_norm_L2                      test_MathUtilities 0)
ADD_TEST( test_norm_L1                      test_MathUtilities 1)
ADD_TEST( test_norm_Linf                    test_MathUtilities 2)
ADD_TEST( test_vec_scale                    test_MathUtilities 3)
ADD_TEST( test_norm_residual_L2             test_MathUtilities 4)
ADD_TEST( test_norm_residual_L1             test_MathUtilities 5)
ADD_TEST( test_norm_residual_Linf           test_MathUtilities 6)

ADD_EXECUTABLE(test_TinyVector              test_TinyVector.cc)
TARGET_LINK_LIBRARIES(test_TinyVector       utilities)
ADD_TEST(test_TinyVector                    test_TinyVector 0)

ADD_EXECUTABLE(test_Random                  test_Random.cc)
TARGET_LINK_LIBRARIES(test_Random           utilities)
ADD_TEST(test_Random                        test_Random 0)

ADD_EXECUTABLE(test_Factory                 test_Factory.cc)
TARGET_LINK_LIBRARIES(test_Factory          utilities)
ADD_TEST(test_Factory                       test_Factory 0)