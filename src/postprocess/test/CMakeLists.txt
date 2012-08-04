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

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_ReactionRates  test_ReactionRates 0)
