//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PulsedExternalSource.hh
 *  @brief  PulsedExternalSource
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_PULSEDEXTERNALSOURCE_HH_
#define detran_PULSEDEXTERNALSOURCE_HH_

#include "TimeDependentExternalSource.hh"

namespace detran
{

/**
 *  @class PulsedExternalSource
 *  @brief External source with Gaussian dependence in time
 *
 *  Given an external source, a Gaussian shape factor in time is
 *  applied.  The Gaussian shape is defined as
 *  @f[
 *      f(t) = e^{- \frac{(t - t_{\text{peak}})^2}{2\sigma^2}}
 *  @f]
 *  where the
 *  @f[
 *      \sigma = \frac{FWHM}{2\sqrt{2 \ln{2}}} \approx 2.23482 \sigma \, .
 *  @f]
 *  The user supplies the time at peak, @f$ t_{\text{peak}} @f$ and
 *  the full width at half maximum, @f$ FWHM @f$.
 *
 *  This shape should be relatively good for modeling pulsed
 *  sources.
 */
class KINETICS_EXPORT PulsedExternalSource: public TimeDependentExternalSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef TimeDependentExternalSource        Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param fixed_source   Fixed source representing peak pulse
   *  @param peak_time      Time of pulse peak
   *  @param fwhm           Full width at half maximum of peak
   *  @param discrete       Flag to indicate treatment as moment or discrete
   */
  PulsedExternalSource(const size_t       number_groups,
                       SP_mesh            mesh,
                       SP_externalsource  fixed_source,
                       const double       peak_time,
                       const double       fwhm,
                       bool               discrete = false);

  /// SP constructor
  static SP_tdsource Create(const size_t       number_groups,
                            SP_mesh            mesh,
                            SP_externalsource  fixed_source,
                            const double       peak_time,
                            const double       fwhm,
                            bool               discrete = false);

  /// Virtual destructor
  virtual ~PulsedExternalSource(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EXTERNAL SOURCES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Get moments source for cell.
   *
   *  Units are n/cc-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   */
  double source(const size_t cell,
                const size_t group);

  /**
   *  @brief Get discrete source for cell and cardinal angle.
   *
   *  Units are n/cc-ster-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   *  @param angle  Cardinal angle index
   */
  double source(const size_t cell,
                const size_t group,
                const size_t angle);

  /**
   *  @brief Set the time of the source
   *  @param time   Time at which source is evaluated
   */
  void set_time(const double time);

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Fixed source
  SP_externalsource d_fixed_source;
  /// Peak time
  double d_peak_time;
  /// Full width at half maximum
  double d_fwhm;
  /// Pulse factor
  double d_factor;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITINS
//---------------------------------------------------------------------------/

#include "PulsedExternalSource.i.hh"

#endif // detran_PULSEDEXTERNALSOURCE_HH_

//---------------------------------------------------------------------------//
//              end of file PulsedExternalSource.hh
//---------------------------------------------------------------------------//
