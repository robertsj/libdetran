#------------------------------------------------------------------------------#
# /src/diffusion/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran diffusion tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#
if(DETRAN_ENABLE_SLEPC)
  ADD_EXECUTABLE(test_DiffusionEigensolver        test_DiffusionEigensolver.cc)
  TARGET_LINK_LIBRARIES(test_DiffusionEigensolver diffusion utilities materials geometry transport)
endif()

if(DETRAN_ENABLE_PETSC)
ADD_EXECUTABLE(test_GainOperator                test_GainOperator.cc)
TARGET_LINK_LIBRARIES(test_GainOperator         diffusion utilities materials geometry)

ADD_EXECUTABLE(test_LossOperator                test_LossOperator.cc)
TARGET_LINK_LIBRARIES(test_LossOperator         diffusion utilities materials geometry)

ADD_EXECUTABLE(test_OneGroupLossOperator        test_OneGroupLossOperator.cc)
TARGET_LINK_LIBRARIES(test_OneGroupLossOperator diffusion utilities materials geometry)
endif()

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

if(DETRAN_ENABLE_SLEPC)
  ADD_TEST( test_DiffusionEigensolver 	test_DiffusionEigensolver 0)
endif()

if(DETRAN_ENABLE_PETSC)
  ADD_TEST( test_GainOperator  	     	test_GainOperator 0)
  ADD_TEST( test_LossOperator  	     	test_LossOperator 0)
  ADD_TEST( test_OneGroupLossOperator  	test_OneGroupLossOperator 0)
endif()