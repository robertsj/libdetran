#------------------------------------------------------------------------------#
# libdetran/src/orthog/test/CMakeLists.txt
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_CLP           test_CLP.cc)
TARGET_LINK_LIBRARIES(test_CLP    orthog)
ADD_TEST(test_CLP   test_CLP      0)

ADD_EXECUTABLE(test_DLP           test_DLP.cc)
TARGET_LINK_LIBRARIES(test_DLP    orthog)
ADD_TEST(test_DLP   test_DLP      0)

ADD_EXECUTABLE(test_DCT           test_DCT.cc)
TARGET_LINK_LIBRARIES(test_DCT    orthog)
ADD_TEST(test_DCT   test_DCT      0)

ADD_EXECUTABLE(test_DCP           test_DCP.cc)
TARGET_LINK_LIBRARIES(test_DCP    orthog)
ADD_TEST(test_DCP   test_DCP      0)

ADD_EXECUTABLE(test_Jacobi01           test_Jacobi01.cc)
TARGET_LINK_LIBRARIES(test_Jacobi01    orthog)
ADD_TEST(test_Jacobi01   test_Jacobi01      0)

ADD_EXECUTABLE(test_ChebyshevU          test_ChebyshevU.cc)
TARGET_LINK_LIBRARIES(test_ChebyshevU   orthog)
ADD_TEST(test_ChebyshevU   test_ChebyshevU      0)

ADD_EXECUTABLE(test_TransformedBasis           test_TransformedBasis.cc)
TARGET_LINK_LIBRARIES(test_TransformedBasis    orthog)
ADD_TEST(test_TransformedBasis                 test_TransformedBasis      0)

ADD_EXECUTABLE(test_BasisComparison           test_BasisComparison.cc)
TARGET_LINK_LIBRARIES(test_BasisComparison    orthog)
ADD_TEST(test_BasisComparisonAnalytic        test_BasisComparison      0)
ADD_TEST(test_BasisComparisonDiscontinuous   test_BasisComparison      1)

ADD_EXECUTABLE(test_BasisIntegration           test_BasisIntegration.cc)
TARGET_LINK_LIBRARIES(test_BasisIntegration    orthog angle)
ADD_TEST(test_BasisIntegration                 test_BasisIntegration      0)

ADD_EXECUTABLE(test_UserBasis           	    test_UserBasis.cc)
TARGET_LINK_LIBRARIES(test_UserBasis    		orthog angle)
ADD_TEST(test_UserBasis                 		test_UserBasis      0)
