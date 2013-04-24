//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Mesh2D.hh
 *  @author Jeremy Roberts
 *  @brief  Mesh2D member definitions.
 */
//---------------------------------------------------------------------------//

#include "Mesh2D.hh"

namespace detran_geometry
{

//---------------------------------------------------------------------------//
Mesh2D::Mesh2D(vec_int xfm, vec_int yfm, 
			   vec_dbl xcme, vec_dbl ycme, 
			   vec_int mat_map)
  : Mesh(2,
         xfm,
         yfm,
         vec_int(1, 1),
         xcme,
         ycme,
         vec_dbl(2, 0.0),
         mat_map)
{
  // Set the second z mesh edge to unity.
  d_zcme[1] = 1.0;
  d_dz[0]   = 1.0;
}

//---------------------------------------------------------------------------//
Mesh2D::Mesh2D(vec_dbl xfme, vec_dbl yfme, vec_int mat_map)
  : Mesh(2,
         xfme,
         yfme,
         vec_dbl(2, 0.0),
         mat_map)
{
  // Set the second z mesh edge to unity.
  d_zcme[1] = 1.0;
  d_dz[0]   = 1.0;
}

} // end namespace detran_geometry

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_IMPLEMENT(detran_geometry::Mesh2D)
#endif

//---------------------------------------------------------------------------//
//              end of Mesh2D.cc
//---------------------------------------------------------------------------//
