//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   angle/MOC.hh
 * \author Joshua J. Jarrell
 * \date   Thu Apr 28 13:34:40 2011
 * \brief  MOC class definition
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: MOC.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef angle_MOC_hh
#define angle_MOC_hh

#include "harness/DBC.hh"
#include "utils/Constants.hh"
#include "Quadrature.hh"

namespace denovo
{

//===========================================================================//
/*!
 * \class MOC
 * \brief 3D MOC quadrature class
 *
 * The MOC quadrature class builds square product quadrature sets for use in a
 * method of characteristics problem. These sets have been developed to ensure
 * for a given region, the rays are consistent across reflecting boundaries in
 * 2D.  It should be noted that these sets have no high order Spherical Harmonic
 * integration capabilities.  
 * 
 * The number of requested azimuthal angles \f$ N_{a,p} \f$ may in general be
 * different than the number of azimuthal angles used in the problem
 * \f$ N_a \f$. This occurs because the rays must reflect correctly on all
 * region boundaries.  The quadrature set determination follows from 
 * Dr. Michael Zika's Ph.D. Dissertation, "Iterative Acceleration For
 * Two-Dimensional Long-Characteristics Transport Problems," Texas A\&M
 * University.
 *
 * We begin with ``provisional'' azimuthal intervals( \f$\phi\f$ is azimuthal
 * angle): 
 * \f[ \Delta\phi_p = \frac{\pi}{2N_{a,p}}\;, \f]
 * and set the angles (\f$m=0:N_{a,r}\f$) to pass through the center of each
 * interval: \f[ \phi_{m,p} = (m+1/2)\Delta\phi_p \; . \f] The provisional
 * angles have a given number of rays per region dimension:
 * \f$ N_{x,m,p} \f$ and \f$ N_{y,m,p} \f$. We can define this using the
 * requested distance between rays, \f$ \delta{a}_{m,p} \f$; the region sizes,
 * \f$\Delta{x}\f$ and \f$\Delta{y}\f$; and the angle in the following manner:
 * \f[ 
    N_{x,m,p} = \frac{\Delta{x}}{\delta{a}_{m,p}}|sin(\phi_{m,p})|
 * \f]
 * and
 * \f[ 
    N_{y,m,p} = \frac{\Delta{y}}{\delta{a}_{m,p}}|cos(\phi_{m,p})|\; .
 * \f]
 * We now use the ceiling function to ensure that the rays align themselves
 * correctly on the region boundaries:
 * \f[ N_{x,m} = \textrm{ceil}(N_{x,m,p}) \f]
 * and
 * \f[ N_{y,m} = \textrm{ceil}(N_{y,m,p}) \;. \f]
 * Using these values results in the following angles:
 * \f[ \phi_{m} = tan^{-1}\left(\frac{\Delta{y}N_{x,m}}
                                     {\Delta{x}N_{y,m}}\right) \; , 
 * \f] 
 * where the arctangent is chosen such that the angle lies in the first
 * quadrant.  The ceiling function can create duplicate angles and we
 * remove the duplicates from the quadrature set.  This results in a ray spacing
 * of 
 * \f[
    \delta{a}_{m} = \frac{\Delta{x}\Delta{y}}
                    { \left(N_{x,m}\Delta{y}\right)^2 + 
                      \left(N_{y,m}\Delta{x}\right)^2  } \;.
 * \f] (Note that the above Equation is incorrect in Zika's Dissertation.)  
 * This results in a distance between rays along the \f$x\f$-axis and 
 * \f$y\f$-axis of 
 * \f[\delta{x}_m=\frac{\delta{a}_{m}}{sin(\phi_{m})} \f] and
 * \f[\delta{y}_m=\frac{\delta{a}_{m}}{cos(\phi_{m})} \f], respectively.
 * The weight of each angle is determined assuming that each angle bisects the
 * interval:
 * \f[
    w_a = \frac{\pi}{2}(\phi_{m+1}-\phi_{m-1}) \;.
 * \f]
 *
 * The polar angles have no requirements to guarantee the reflecting boundaries.
 * The polar intervals are defined \f[ \Delta{\Theta}_{p} =\frac{\pi}{2N_p}\;,\f]
 * where \f$N_p\f$ is the number of polar levels and \f$\Theta_{p}\f$ is the
 * polar angle.  The polar directions are chosen to preserve the mean secant:
 * \f[
    sec(\Theta_p)=\frac{
         \int_{\Theta_{p-1/2}}^{\Theta_{p+1/2}}sec(\Theta)cos^2(\Theta)d\Theta}
        {\int_{\Theta_{p-1/2}}^{\Theta_{p+1/2}}cos^2(\Theta)d\Theta} \; ,
 * \f]
 * where the edges of the polar intervals are 
 * \f$ \Theta_{p-1/2} = (p-1)\Delta{\Theta}_p \f$ and 
 * \f$ \Theta_{p+1/2} = p\Delta{\Theta}_p \f$.  By preserving the mean secant, we
 * can accurately compute the first-flight collision probability. The polar
 * angles are
 * \f[
    \Theta_p = cos^{-1}\left( \frac{ \frac{\pi}{4N_p} +
        \frac{1}{4}\left(sin(2\Theta_{p+1/2})-sin(2\Theta_{p-1/2})\right)}
        {sin(\Theta_{p+1/2})-sin(\Theta_{p-1/2})}\right)
 * \f]
 * and their associated weights are:
 * \f[
    w_p = \int_{\Theta_{p-1/2}}^{\Theta_{p+1/2}}d\Theta{cos(\Theta)}
        = sin(\Theta_{p+1/2})-sin(\Theta_{p-1/2}) \;.
 * \f]
 *
 * Finally, product quadrature rule establishes that direction
 * for each ordinate is
 * \f[ \mu  = cos(\phi)cos(\Theta) \; , \\
       \eta = sin(\phi)cos(\Theta) \; , \\
       \xi  = sin(\Theta) \; , 
 * \f]
 * and the weight is the product of the polar weight and the azimuthal weight, 
 * \f$ w = w_a{w_p} \f$.  We then reflect the first octant sets into the other
 * seven octants to completely define the quadrature set.
 *
 * We note that this is only one option for creating a method of characteristics
 * quadrature set.  Other sets with higher order integration capabilities have
 * been derived. 
 */
/*! 
 * \example angle/test/tstMOC.cc
 *
 * Test of Method of Characteristics Quadrature Sets
 */
//===========================================================================//

class MOC : public Quadrature
{
  private:
    //! Base class typdef
    typedef Quadrature Base;

  public:

    // Constructor for square set
    MOC( int num_levels, double delta_x, double delta_y, double ray_spacing,
         double norm=4*constants::pi);
    
    //>>> INHERITED INTERFACE
    
    // Display the quadrature.
    void display() const;

    //! Label.
    std_string label() const { return "moc3D"; }

    //! Get the weights.
    const Vec_Dbl& weight() const { return b_wt; }

    //! Return the weight for a given ordinate.
    double weight(int m) const { Require (m < b_wt.size()); return b_wt[m]; }
};

} // end namespace denovo

#endif // angle_MOC_hh

//---------------------------------------------------------------------------//
//              end of angle/MOC.hh
//---------------------------------------------------------------------------//
