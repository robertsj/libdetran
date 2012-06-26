#------------------------------------------------------------------------------#
# libdetran/src/moc/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran moc tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

# Azimuthal Quadratures
ADD_EXECUTABLE(test_Collocated 	test_Collocated.cc)
ADD_EXECUTABLE(test_Uniform 	test_Uniform.cc)
# Polar Quadratures
ADD_EXECUTABLE(test_TabuchiYamamoto test_TabuchiYamamoto.cc)
# Other
ADD_EXECUTABLE(test_Point   test_Point.cc)
ADD_EXECUTABLE(test_Segment test_Segment.cc)
ADD_EXECUTABLE(test_Track   test_Track.cc)
ADD_EXECUTABLE(test_Tracker test_Tracker.cc)

#
TARGET_LINK_LIBRARIES(test_Collocated      moc angle utilities)
TARGET_LINK_LIBRARIES(test_Uniform         moc angle utilities)
#
TARGET_LINK_LIBRARIES(test_TabuchiYamamoto moc angle utilities)
#
TARGET_LINK_LIBRARIES(test_Point           moc angle utilities)
TARGET_LINK_LIBRARIES(test_Segment         moc angle utilities)
TARGET_LINK_LIBRARIES(test_Track           moc angle utilities)
TARGET_LINK_LIBRARIES(test_Tracker         moc angle utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_Collocated      test_Collocated 0)
ADD_TEST(test_Uniform         test_Uniform    0)
#
ADD_TEST(test_TabuchiYamamoto test_TabuchiYamamoto 0)
#
ADD_TEST(test_Point           test_Point   0)
ADD_TEST(test_Segment         test_Segment 0)
ADD_TEST(test_Track           test_Track   0)
ADD_TEST(test_Tracker_2x2     test_Tracker 0)
ADD_TEST(test_Tracker_3x3     test_Tracker 1)


