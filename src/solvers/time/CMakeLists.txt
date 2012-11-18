# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TIME_SRC
  ${SRC_DIR}/TimeStepper.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/solver/time .hh)

