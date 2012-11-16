//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeDependentExternalSource.hh
 *  @brief  TimeDependentExternalSource
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMEDEPENDENTEXTERNALSOURCE_HH_
#define detran_TIMEDEPENDENTEXTERNALSOURCE_HH_

#include "external_source/ExternalSource.hh"

namespace detran
{

/**
 *  @class TimeDependentExternalSource
 *  @brief Base class for time-dependent external sources
 *
 *  The only addition beyond the base external source class is to
 *  specify a time at which the source will be evaluated during the
 *  transport solve for the current step.
 */
class TimeDependentExternalSource: public detran_external_source::ExternalSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_external_source::ExternalSource        Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param quadrature     Pointer to angular quadrature
   */
  TimeDependentExternalSource(size_t        number_groups,
                              SP_mesh       mesh,
                              SP_quadrature quadrature)
    : Base(number_groups, mesh, quadrature)
    , d_time(0.0)
  {
    /* ... */
  }

  /// Virtual destructor
  virtual ~TimeDependentExternalSource(){}

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
  virtual double source(const size_t cell,
                        const size_t group) = 0;

  /**
   *  @brief Get discrete source for cell and cardinal angle.
   *
   *  Units are n/cc-ster-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   *  @param angle  Cardinal angle index
   */
  virtual double source(const size_t cell,
                        const size_t group,
                        const size_t angle) = 0;

  /**
   *  @brief Set the time of the source
   *  @param time   Time at which source is evaluated
   */
  virtual void set_time(const double time) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Current time [s]
  double d_time;

};

} // end namespace detran

#endif // detran_TIMEDEPENDENTEXTERNALSOURCE_HH_

//---------------------------------------------------------------------------//
//              end of file TimeDependentExternalSource.hh
//---------------------------------------------------------------------------//
