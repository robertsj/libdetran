#------------------------------------------------------------------------------#
# libdetran/src/moc/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(MOC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(MOC_SRC
	QuadratureMOC.cc
    TabuchiYamamoto.cc
    TrackDB.cc
    Tracker.cc
    Uniform.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(moc
            SHARED
            ${MOC_SRC} 
)
TARGET_LINK_LIBRARIES( moc
					   geometry
					   materials
                       utilities
                       ${GPERFTOOLS_LIBRARIES}                        
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib moc)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${MOC_DIR}/test)
    add_subdirectory(${MOC_DIR}/test)
endif()
