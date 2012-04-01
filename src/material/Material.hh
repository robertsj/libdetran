//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Material.hh
 * \author Jeremy Roberts
 * \brief  Material class definition.
 */
//---------------------------------------------------------------------------//


#ifndef MATERIAL_HH_
#define MATERIAL_HH_

// Other detran headers
#include "Definitions.hh"
#include "DBC.hh"
#include "SP.hh"
#include "Warning.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class Material
 * \brief Simple cross section container.
 */
//===========================================================================//
class Material : public detran_utils::Object
{

public:

  /// \name Useful typedefs when using Material
  /// \{
  typedef typename detran_utils::SP<Material>   SP_material;
  typedef detran_utils::vec_dbl                 vec_dbl;
  typedef detran_utils::vec_int                 vec_int;
  typedef detran_utils::vec2_dbl                vec2_dbl;
  typedef detran_utils::vec2_int                vec2_int;
  typedef detran_utils::vec3_dbl                vec3_dbl;
  /// \}

  /*!
   *  \brief Constructor.
   *
   *  \param    number_groups       Number of energy groups.
   *  \param    number_materials    Number of materials.
   *  \param    downscatter         Switch on to use only downscatter.
   */
  Material(int  number_groups,
           int  number_materials,
           bool downscatter = false);

  /*!
   *  \brief SP Constructor.
   *
   *  \param    number_groups       Number of energy groups.
   *  \param    number_materials    Number of materials.
   *  \param    downscatter         Switch on to use only downscatter.
   *  \return                       Smart pointer to Material object.
   */
  static SP_material Create(int number_groups,
                            int number_materials,
                            bool downscatter)
  {
    SP_material p;
    p = new Material(number_groups, number_materials, downscatter);
    return p;
  }

  //--------------------------------------------------------------------------//
  // Setters
  //--------------------------------------------------------------------------//

  void set_sigma_t(int m, int g, double v);

  void set_sigma_a(int m, int g, double v);

  void set_nu_sigma_f(int m, int g, double v);

  void set_chi(int m, int g, double v);

  void set_sigma_s(int m, int g, int gp, double v);

  void set_diff_coef(int m, int g, double v);

  // Vectorized setters

  void set_sigma_t(int m, vec_dbl &v);

  void set_nu_sigma_f(int m, vec_dbl &v);

  void set_chi(int m, vec_dbl &v);

  void set_sigma_s(int m, int g, vec_dbl &v);

  void set_diff_coef(int m, vec_dbl &v);

  //------------------------------------------------------------------------//
  // Getters
  //------------------------------------------------------------------------//

  inline double sigma_t(int m, int g) const;

  inline double sigma_a(int m, int g) const;

  inline double nu_sigma_f(int m, int g) const;

  inline double chi(int m, int g) const;

  inline double sigma_s(int m, int g, int gp) const;

  inline double diff_coef(int m, int g) const;

  int number_groups()
  {
    return d_number_groups;
  }

  int number_materials()
  {
    return d_number_materials;
  }

  /*!
   *  \brief Lower scatter group bound.
   *
   *  This is the *lowest* index (highest energy) \f$ g' \f$
   *  that leads to downscatter for a given outgoing group \f$ g \f$.
   */
  int lower(int g);

  /*!
   *  \brief Upper scatter group bound.
   *
   *  This is the *highest* index (lowest energy) \f$ g' \f$
   *  that upscatters into the outgoing group \f$ g \f$.
   */
  int upper(int g);

  /// Do we do only downscatter?
  bool downscatter()
  {
    return d_downscatter;
  }

  /// Index below which upscatter doesn't occur for any material.
  int upscatter_cutoff()
  {
    return d_upscatter_cutoff;
  }

  /// Computes scattering bounds and absorption cross section.
  void finalize();

  /// Pretty print the material database.
  void display();

  /// Incomplete implementation of DBC function.
  bool is_valid() const
  {
    Ensure(d_number_groups >= 0);
    Ensure(d_finalized);
  };

private:

  /// Number of groups
  int d_number_groups;

  /// Number of materials
  int d_number_materials;

  /// Downscatter switch (when true, upscatter ignored)
  bool d_downscatter;

  /// Total cross section [material, group]
  vec2_dbl d_sigma_t;

  /// Absorption cross section [material, group]
  vec2_dbl d_sigma_a;

  /// Fission [material, group]
  vec2_dbl d_nu_sigma_f;

  /// Fission spectrum [material, group]
  vec2_dbl d_chi;

  /// Scatter [material, group<-, group']
  vec3_dbl d_sigma_s;

  /// Diffusion coefficient [material, group]
  vec2_dbl d_diff_coef;

  /// Scatter bounds applied to all materials [group, 2]
  vec2_int d_scatter_bounds;

  /*!
   *  Upscatter cutoff.  Only groups equal to or above this cutoff are
   *  subject to upscatter iterations.
   */
  int d_upscatter_cutoff;

  /// Are we ready to be used?
  bool d_finalized;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Material.i.hh"

#endif /* MATERIAL_HH_ */
