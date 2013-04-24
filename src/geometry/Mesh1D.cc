//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Mesh1D.hh
 *  @author Jeremy Roberts
 *  @brief  Mesh1D member definitions.
 */
//---------------------------------------------------------------------------//

#include "Mesh1D.hh"

namespace detran_geometry
{

//---------------------------------------------------------------------------//
Mesh1D::Mesh1D(vec_int xfm, vec_dbl xcme, vec_int mat_map)
  : Mesh(1,
         xfm,
         vec_int(1, 1),
         vec_int(1, 1),
         xcme,
         vec_dbl(2, 0.0),
         vec_dbl(2, 0.0),
         mat_map)
{
  // Set the second y mesh edge to unity.
  d_ycme[1] = 1.0;
  d_dy[0]   = 1.0;
  // Set the second z mesh edge to unity.
  d_zcme[1] = 1.0;
  d_dz[0]   = 1.0;
}

//---------------------------------------------------------------------------//
Mesh1D::Mesh1D(vec_dbl xfme, vec_int mat_map)
  : Mesh(1,
         xfme,
         vec_dbl(2, 0.0),
         vec_dbl(2, 0.0),
         mat_map)
{
  // Set the second y mesh edge to unity.
  d_ycme[1] = 1.0;
  d_dy[0]   = 1.0;
  // Set the second z mesh edge to unity.
  d_zcme[1] = 1.0;
  d_dz[0]   = 1.0;
}

} // end namespace detran_geometry

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_IMPLEMENT(detran_geometry::Mesh1D)
#endif

//---------------------------------------------------------------------------//
//              end of Mesh1D.cc
//---------------------------------------------------------------------------//
