//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearMaterial.hh
 *  @brief  LinearMaterial
 *  @author Jeremy Roberts
 *  @date   Nov 14, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_LINEARMATERIAL_HH_
#define detran_LINEARMATERIAL_HH_

#include "TimeDependentMaterial.hh"

namespace detran
{

/**
 *  @class LinearMaterial
 *  @brief Material defined via discrete materials to be linearly
 *         interpolated
 *
 *  The user specified material databases for a sequence of monotonically
 *  increasing times,
 *  @f[
 *      t_0, t_1, \cdots, t_N \quad t_0 \ge 0
 *  @f]
 *  If the initial time is not zero, then for times
 *  times @f$ t \le t_0 @f$, the first material is
 *  returned.  Likewise, if $t \ge t_N$, the last
 *  material is returned.
 */
class KINETICS_EXPORT LinearMaterial: public TimeDependentMaterial
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef TimeDependentMaterial             Base;
  typedef KineticsMaterial::SP_material     SP_kineticsmaterial;
  typedef std::vector<SP_kineticsmaterial>  vec_material;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param times                      Times at which materials are defined
   *  @param materials                  Materials defined at discrete times
   *  @note The first time *must* be zero
   */
  LinearMaterial(const vec_dbl      &times,
                 const vec_material &materials,
                 std::string        name = "LinearMaterial");

  /// SP constructor
  static SP_material Create(const vec_dbl      &times,
                            const vec_material &materials,
                            std::string        name = "LinearMaterials");

  /// Virtual destructor
  virtual ~LinearMaterial(){};

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Times at which materials are defined
  vec_dbl d_times;
  /// Number of times
  size_t d_number_times;
  /// Materials at discrete times
  vec_material d_materials;
  /// Temporary material container for interpolated values
  SP_kineticsmaterial d_temporary;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL TIME DEPENDENT MATERIALS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Linearly interpolate materials
  void update_impl();

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /**
   *  @brief Update the material via interpolation
   *  @param index_A  Time index of first material
   *  @param index_B  Time index of second material
   *  @param A        Coefficient of first material
   *  @param B        Coefficient of second material
   */
  void update_material(const size_t index_A,
                       const size_t index_B,
                       const double A,
                       const double B);

};

} // end namespace detran

#endif // detran_LINEARMATERIAL_HH_ 

//---------------------------------------------------------------------------//
//              end of file LinearMaterial.hh
//---------------------------------------------------------------------------//
