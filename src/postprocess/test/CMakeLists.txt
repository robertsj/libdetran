#------------------------------------------------------------------------------#
# libdetran/src/postprocess/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran postprocess tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_ReactionRates test_ReactionRates.cc)
TARGET_LINK_LIBRARIES(test_ReactionRates postprocess 
                                      	 geometry 
                                         materials 
                                         transport 
                                         utilities)

if(DETRAN_ENABLE_SILO)
ADD_EXECUTABLE(test_SiloOutput test_SiloOutput.cc)
TARGET_LINK_LIBRARIES(test_SiloOutput postprocess 
                                      geometry 
                                      materials 
                                      transport 
                                      utilities 
                                      ${Silo_LIBRARIES})
endif()

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_ReactionRates  test_ReactionRates 0)

if(DETRAN_ENABLE_SILO)
ADD_TEST( test_SiloOutput  test_SiloOutput 0)
endif()