//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearExternalSource.hh
 *  @brief  LinearExternalSource
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_LINEAREXTERNALSOURCE_HH_
#define detran_LINEAREXTERNALSOURCE_HH_

#include "TimeDependentExternalSource.hh"

namespace detran
{

/**
 *  @class LinearExternalSource
 *  @brief Base class for time-dependent external sources
 *
 *  The only addition beyond the base external source class is to
 *  specify a time at which the source will be evaluated during the
 *  transport solve for the current step.
 */
class KINETICS_EXPORT LinearExternalSource: public TimeDependentExternalSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef TimeDependentExternalSource        Base;
  typedef std::vector<SP_externalsource>     vec_source;
  typedef detran_utilities::vec_dbl          vec_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param times          Times at which sources are defined
   *  @param sources        Sources at each time
   *  @param discrete       Flag to indicate treatment as moment or discrete
   */
  LinearExternalSource(const size_t   number_groups,
                       SP_mesh        mesh,
                       vec_dbl        times,
                       vec_source     sources,
                       bool           discrete = false);

  /// SP constructor
  static SP_tdsource Create(const size_t   number_groups,
                            SP_mesh        mesh,
                            vec_dbl        times,
                            vec_source     sources,
                            bool           discrete = false);

  /// Virtual destructor
  virtual ~LinearExternalSource(){}

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

  /// Discrete times at which sources are defined
  vec_dbl d_times;
  /// Number of times
  double d_number_times;
  /// Sources at each time
  vec_source d_sources;
  /// Index of first source
  size_t d_ia;
  /// Interpolation factor for first source
  double d_fa;
  /// Index of second source
  size_t d_ib;
  /// Interpolation factor for second source
  double d_fb;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITINS
//---------------------------------------------------------------------------/

#include "LinearExternalSource.i.hh"

#endif // detran_LINEAREXTERNALSOURCE_HH_ 

//---------------------------------------------------------------------------//
//              end of file LinearExternalSource.hh
//---------------------------------------------------------------------------//
