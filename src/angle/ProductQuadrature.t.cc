//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.t.cc
 *  @brief ProductQuadrature explicit instantiations
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include "angle/ProductQuadrature.t.hh"

namespace detran_angle
{

// product
//
template class ProductQuadrature<AzimuthalGL, PolarGL>;
template class ProductQuadrature<AzimuthalGL, PolarDGL>;
template class ProductQuadrature<AzimuthalGL, PolarU>;
template class ProductQuadrature<AzimuthalGL, PolarTY>;
//
template class ProductQuadrature<AzimuthalDGL, PolarGL>;
template class ProductQuadrature<AzimuthalDGL, PolarDGL>;
template class ProductQuadrature<AzimuthalDGL, PolarU>;
template class ProductQuadrature<AzimuthalDGL, PolarTY>;
//
template class ProductQuadrature<AzimuthalU, PolarGL>;
template class ProductQuadrature<AzimuthalU, PolarDGL>;
template class ProductQuadrature<AzimuthalU, PolarU>;
template class ProductQuadrature<AzimuthalU, PolarTY>;
//
template class ProductQuadrature<AzimuthalS, PolarGL>;
template class ProductQuadrature<AzimuthalS, PolarDGL>;
template class ProductQuadrature<AzimuthalS, PolarU>;
template class ProductQuadrature<AzimuthalS, PolarTY>;
//
template class ProductQuadrature<AzimuthalASQR, PolarASDR>;


} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.t.cc
//----------------------------------------------------------------------------//
