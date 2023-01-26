//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeIndependentMaterial.hh
 *  @author Rabab Elzohery
 *  @date   March 28, 2020
 *  @brief  TimeIndependentMaterial class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMEINDEPENDENTMATERIAL_HH_
#define detran_TIMEINDEPENDENTMATERIAL_HH_

#include "TimeDependentMaterial.hh"
#include "utilities/Definitions.hh"
#include <memory>

namespace detran
{

/**
 *  @class TimeIndependentMaterial
 *  @brief Base class for time-Independent materials
 */
class KINETICS_EXPORT TimeIndependentMaterial: public TimeDependentMaterial
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef TimeDependentMaterial                             Base;
  typedef Base::SP_material                                 SP_base;
  typedef std::shared_ptr<TimeIndependentMaterial>	SP_material;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param d_km          Referene of KineticsMaterial
   */
  TimeIndependentMaterial(KineticsMaterial::SP_material km);

  static SP_material Create(KineticsMaterial::SP_material km);

  KineticsMaterial ::SP_material d_km;

  void update_impl();

  /// Virtual destructor
  virtual ~TimeIndependentMaterial(){}

};

KINETICS_TEMPLATE_EXPORT(detran_utilities::SP<TimeIndependentMaterial>)

} // end namespace detran

#endif /* detran_TIMEDEPENDENTMATERIAL_HH_ */
