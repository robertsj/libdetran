#------------------------------------------------------------------------------#
# src/angle/CMakeLists.txt
#------------------------------------------------------------------------------#

set(LIB_NAME "angle")

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

set(SRC
    AbuShumaysDoubleRange.cc
    AbuShumaysQuadrupleRange.cc
    #Collocated.cc
    LevelSymmetric.cc
    MomentToDiscrete.cc
    Quadrature.cc
    QuadratureFactory.cc
    #QuadratureMOC.cc
    #QuadrupleRange.cc
    SphericalHarmonics.cc
    TabuchiYamamoto.cc
    #Uniform.cc
    UniformEqual.cc
    #ChebyshevLegendre.cc
    #ChebyshevDPN.cc
    MomentIndexer.cc
    #LegendreDTN.cc
    PolarQuadrature.t.cc
    AzimuthalQuadrature.t.cc
    ProductQuadrature.cc
    ProductQuadratureAP.t.cc
)

#------------------------------------------------------------------------------#
# LINKED LIBRARIES
#------------------------------------------------------------------------------#

set(LINKED_LIBS
    utilities
)

#------------------------------------------------------------------------------#
# LIBRARY
#------------------------------------------------------------------------------#

setup_detran_library(${LIB_NAME} "${SRC}" "${LINKED_LIBS}")
