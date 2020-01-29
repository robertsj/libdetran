//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  AzimuthalQuadrature.t.cc
 *  @brief AzimuthalQuadrature explicit instantiations
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "angle/AzimuthalQuadrature.t.hh"

namespace detran_angle
{

template class AzimuthalQuadrature<BaseGL>;
template class AzimuthalQuadrature<BaseDGL>;
template class AzimuthalQuadrature<BaseUniform>;
template class AzimuthalQuadrature<BaseSimpson>;
template class AzimuthalQuadrature<AbuShumaysQuadrupleRange>;

} // end namespace detran_angle



