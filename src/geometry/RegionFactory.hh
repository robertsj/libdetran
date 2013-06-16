//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  RegionFactory.hh
 *  @brief RegionFactory struct defition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_REGIONFACTORY_HH_
#define detran_geometry_REGIONFACTORY_HH_

#include "Region.hh"
#include "QuadraticSurface.hh"
#include "QuadraticSurfaceFactory.hh"
#include "PinCell.hh"
#include "Assembly.hh"

namespace detran_geometry
{

/**
 *  @struct RegionFactory
 *  @brief  Convenience functions for creating common regions
 */
struct RegionFactory
{

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//



  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef PinCell::SP_pincell               SP_pincell;
  typedef Assembly::SP_assembly             SP_assembly;
  typedef Region::SP_region                 SP_region;
  typedef Region::vec_region                vec_region;
  typedef QuadraticSurface                  QS;
  typedef Surface::SP_surface               SP_surface;
  typedef Surface::vec_surface              vec_surface;
  typedef Point                             c_Point;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec_int         vec_int;

  typedef QuadraticSurfaceFactory           QF;


  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Create regions from a pin cell
   *  @param   pin      Pin cell
   *  @return           Region representation of pin cell
   */
  static vec_region CreatePinCell(SP_pincell pin);

  /**
   *  @brief Create a region from an assembly
   *  @param   pin      Pin cell
   *  @return           Region representation of pin cell
   */
  static vec_region CreateAssembly(SP_assembly  assembly);

  /// Translate a region to a new origin
  static SP_region CreateTranslatedRegion(SP_region r, const Point &origin);

};

} // end namespace detran_geometry

#endif /* detran_geometry_REGIONFACTORY_HH_ */

//----------------------------------------------------------------------------//
//              end of RegionFactory.hh
//----------------------------------------------------------------------------//
