//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Collocated.hh
 * \brief  Collocated 
 * \author Jeremy Roberts
 * \date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef COLLOCATED_HH_
#define COLLOCATED_HH_

#include "QuadratureMOC.hh"

namespace detran_angle
{

/*!
 *  \class Collocated
 *  \brief A MOC quadrature with a finite set of origins along a side.
 *
 * This quadrature has a number of unique features.  It satisfies cyclic
 * tracking requirements in a rectangular region.  It also uses a finite
 * set of points along a surface as track entrance (and exit) points. In
 * this first implementation, we consider only square regions.
 *
 * The user sets the number of spatial points along a side of a unit
 * cell.  Currently,
 * this must be a power of 3, and only 9 and 27 are supported.
 * Additionally, a multiple may be assigned so that the base points
 * are repeated on adjacent unit cells.  This lets us define the
 * spacing in terms of a pin cell, which is then repeated on
 * the assembly level.
 *
 * The angles (between \f$ \pi/4 \f$ and \f$ \pi/2 \f$) are uniquely defined
 * by
 * \f[
 *      \tan{\phi_i} = 3^i \, , \,\,\, i = 0, \, \cdots n \, ,
 * \f]
 * where
 * \f[
 *      n = \mathrm{floor}(\log_3(N)+1) \,
 * \f]
 * and  \em N is the number of spatial points. Note that n/N is
 * maximized when \em N  is a power of three.  Using 9 spatial points
 * surface yields 5 angles per quadrant.  27 points yields 7 angles per
 * quadrant.  The other angles are defined by symmetry about \f$ \pi/4 \f$.
 *
 * For weights, the user can use an arc length approximation or a
 * quadrature rule defined for the points that may provide better
 * accuracy.
 *
 */
class ANGLE_EXPORT Collocated: public QuadratureMOC
{

public:

  /*!
   *  \brief Constructor
   *  \param dim                  Problem dimension (only 2 supported for now)
   *  \param num_azimuths_octant  Number of azimuths per octant
   *  \param num_space            Number of tracks per azimuth
   *  \param num_polar            Number of polar angles in half space
   *  \param polar                Polar quadrature string identifier
   */
  Collocated(size_t dim,
             size_t num_azimuths_octant,
             size_t multiplier,
             size_t num_polar,
             std::string polar);

  /// SP Constructor.
  static detran_utilities::SP<Quadrature>
    Create(size_t dim, size_t num_azimuths_octant, size_t multiplier,
           size_t num_polar, std::string polar)
  {
    SP_quadrature p(new Collocated(dim, num_azimuths_octant,
                                   multiplier, num_polar, polar));
    return p;
  }

  ~Collocated(){};

private:


};

} // end namespace detran_angle

#endif // COLLOCATED_HH_ 

//---------------------------------------------------------------------------//
//              end of file Collocated.hh
//---------------------------------------------------------------------------//
