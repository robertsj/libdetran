//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SphericalHarmonics.hh
 *  @author Jeremy Roberts
 *  @date   Jun 29, 2011
 *  @brief  SphericalHarmonics class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_SPHERICAL_HARMONICS_HH_
#define detran_angle_SPHERICAL_HARMONICS_HH_

#include "angle/angle_export.hh"
#include "utilities/Definitions.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
/**
 *  @class SphericalHarmonics
 *  @brief Spherical harmonics generation for anisotropic scattering.
 *
 *  In 2-d and 3-d problems, an expansion of the angular flux in spherical
 *  harmonics is required to treat anisotropic scattering.  In 1-d problems,
 *  the angular flux is expanded in Legendre polynomials, which
 *  can be considered to be a subset of the spherical harmonics.
 *
 *  We follow the presentation of Hebert.
 *
 *  A function of a directional cosine \f$\xi\f$, \f$ f(\xi) \f$, can
 *  be represented as in \f$L\f$-th order Legendre expansion:
 *  \f[
 *      f( \xi) \approx \sum^L_{l=0} \frac{2l+1}{2}P_l(\xi)f_l \, ,
 *  \f]
 *  where the Legendre coefficients are defined
 *  \f[
        f_l = \int^{1}_{-1} d\xi P_l(\xi) f(\xi) \, .
 *  \f]
 *  In the standard 1-d formulation, in which the azimuthal dependence is
 *  removed by integration, the zeroth Legendre moment angular
 *  flux \f$ \psi(z,\xi) \f$ is
 *  \f[
 *      \phi_0(z) = \int^{1}_{-1} d\xi  \psi(z,\xi) \, ,
 *  \f]
 *  which we recognize as the scalar flux.  The scattering source is
 *  defined in 1D as
 *  \f[
 *      Q(z,\xi) = \int^{2\pi}_{0} d\phi' \int^{1}_{-1}d\xi'
 *                      \Sigma_s(z,\mu_0)\psi(z,\xi') \\
 *  \f]
 * and can likewise be expanded (skipping several important
 * steps! see e.g. the 22106 notes) as
 *  @f[
       Q(z,\xi) \approx \sum^{L}_{l=0}
                      \frac{2l+1}{2}\Sigma_{sl}(z)P_l(\xi)\phi_l(z) \, ,
 *  @f]
 * where the Legendre moments \f$ \Sigma_{sl} \f$ are provided in
 * a cross-section library.  Note, for consistency with the
 * multidimensional coordinates, we place 1D problems along
 * the \f$z\f$-axis.  For 1D problems using an \$L\$-th order
 * expansion (for scattering), there are \f$L+1\f$ flux moments
 * per node (and possibly more than one node per mesh).
 *
 * For 2D and 3D problems, the azimuthal dependence requires use of
 * the spherical harmonics for expansions.
 *
 * The spherical harmonics \f$Y^m_l(\Omega) = Y^m_l(\xi,\varphi)\f$
 * (in the source, denoted \f$Y_{lm}\f$ ) are defined
 *  @f[
       Y^m_l(\xi,\varphi) = \sqrt{ (2-\delta_{m,0})
                                      \frac{(l-|m|)!}{(l+|m|)!} }
                               P^{|m|}_l (\xi) \mathcal{T}_m (\varphi ) \ .
 *  @f]
 * where the associated Legendre polynomials are defined in terms of
 * Legendre polynomials as
 *  @f[
       P^m_l(\xi) = (1-\xi^2)^{m/2} \frac{d^m}{d\xi^m}
                     P_l(\xi) \, , \,\,\,\,\, m \geq 0 \, ,
 *  @f]
 * and the trigonometric functions are
 *  @f[
       \mathcal{T}_m (\varphi ) \left\{
         \begin{array}{l l}
            \cos(m\varphi)      & \quad \text{if $m \geq 0$ }\\
            \sin(|m|\varphi)    & \quad \text{otherwise} \, . \\
         \end{array} \right.
 *  @f]
 * The associated Legendre polynomials are as given use the so-called Ferrer
 * definition, omitting the typically standard factor of \f$ (-1)^m \f$.
 * Hebert notes this representation is helpful since
 *  @f[
         \mathbf{\Omega} =  \begin{pmatrix}
        \mu      \\
        \eta     \\
        \xi
     \end{pmatrix} =  \begin{pmatrix}
        \sin(\theta)\cos(\varphi)      \\
        \sin(\theta)\sin(\varphi)      \\
        \cos(\theta)
     \end{pmatrix} = \begin{pmatrix}
        Y^{1}_{1}(\xi,\varphi)   \\
        Y^{-1}_{1}(\xi,\varphi)  \\
        Y^{0}_{1}(\xi,\varphi)
       \end{pmatrix} \, .
 *  @f]
 * The spherical harmonics through \f$L=3\f$ in terms of the directional
 * cosines are
 *  @f[
 *     Y^{0}_{0} = 1 \, ,
 *  @f]
 *  @f[
 *     Y^{-1}_{1}=\eta \, , \quad
 *     Y^{0}_{1}=\xi   \, , \quad
 *     Y^{1}_{1}=\mu
 *  @f]
 *  @f[
 *     Y^{-2}_{2} = \sqrt{3}\mu\eta          \, , \quad
 *     Y^{-1}_{2}=\sqrt{3}\xi\eta            \, , \quad
 *     Y^{0}_{2}=\frac{1}{2}(3\xi^2-1)       \, , \quad
 *     Y^{1}_{2}=\sqrt{3}\xi\mu              \, , \quad
 *     Y^{2}_{2}=\frac{\sqrt{3}}{2}(\mu^2-\eta^2) \, .
 *  @f]
 *  @f[
 *     Y^{-3}_{3} = \sqrt{\frac{5}{8}}\eta(3\mu^2-\eta^2)      \, , \quad
 *     Y^{-2}_{3} = \sqrt{15}\xi\mu\eta                        \, , \quad
 *     Y^{-1}_{3} = \sqrt{\frac{3}{8}}\eta(5\xi^2-1)           \, , \quad
 *     Y^{0}_{3}  = \frac{1}{2}(5\xi^3-3\xi)                   \, , \quad
 *     Y^{1}_{3}  = \sqrt{\frac{3}{8}}\mu(5\xi^2-1)            \, , \quad
 *     Y^{2}_{3}  = \sqrt{\frac{15}{4}}\xi(\mu^2-\eta^2)       \, , \quad
 *     Y^{3}_{3}  = \sqrt{\frac{5}{8}}\mu(\mu^2-3\eta^2)       \, .
 *  @f]
 * The trigonometric functions, associated Legendre polynomials, and
 * spherical harmonics satisfy the following orthogonality conditions:
 *  @f[
        \int^{\pi}_{-\pi} d\varphi \mathcal{T}_n (\varphi )
           \mathcal{T}_{n'}(\varphi ) = \pi(1+\delta_{n,0})\delta_{n,n'} \, ,
 *  @f]
 *  @f[
        \int^{1}_{-1} d\mu P^m_l (\mu)
          P^{m'}_{l'} = \frac{2(l+m)!}{(2l+1)(l-m)!}\delta_{l,l'} \, ,
 *  @f]
 * and
 *  @f[
         \int_{4\pi} d^2\Omega Y^m_l(\mathbf{\Omega})Y^m_{l'}(\mathbf{\Omega})
             = \int^{2\pi}_{0} d\phi' \int^{1}_{-1}d\mu'
               Y^m_l(\theta,\varphi) Y^{m'}_{l'}(\theta,\varphi)
             = \frac{4\pi}{2l+1} \delta_{l,l'}\delta_{m,m'} \, .
 *  @f]
 * Then, angular flux is expanded following
 *  @f[
 *     \psi(\mathbf{r},\mathbf{\Omega}) \approx
 *       \sum^{L}_{l=0} \frac{2l+1}{4\pi}
 *       \sum^{l}_{m=-l} Y^m_l (\mathbf{\Omega}) \phi^m_l(\mathbf{r}) \, ,
 *  @f]
 * where
 *  @f[
       \phi^m_l(\mathbf{r}) = \int_{4\pi} d^2\Omega Y^m_l(\mathbf{\Omega})
                              \psi(\mathbf{r},\mathbf{\Omega}) \, .
 *  @f]
 * The additional \f$2\pi\f$ in the expansion comes from the isotropy in
 * the azimuthal angle, multiplied away in 1D (a 1D angular flux is really
 * in units of 1/rad, not 1/steradian).
 *
 * Again, we note that the zeroth order moment is the angular flux.  For
 * an \f$L\f$th order expansion, there are \f$(L+1)^2\f$ moments.
 *
 *  @note
 * For 2-D, the moments \f$\phi^0_{l>0}\f$
 * vanish identically.  This is because a 2D problem is defined such
 * that there is no variation in the polar (\f$z\f$) direction.  The moments
 *  @f$\phi^0_{l>0}\f$ represent net changes in the \f$z\f$ direction, and
 * so should be eliminated for efficiency.  For 2D problems, the number of
 * moments is therefore \f$(L+1)^2-L\f$. This is handled in
 *  @ref Moment_to_Discrete.
 *
 * It should be noted that the for \f$m=0 \f$, the spherical harmonics as
 * defined reduce to the Legendre polynomials in \f$\xi\f$.  Hence, this
 * class can be used for 1D expansions.
 *
 * Currently, only the spherical harmonics through \f$ L=3 \f$ are
 * implemented, though adding capability for arbitrary orders should be
 * straightforward if external libraries are used.
 *
 *  @note
 * Legendre expansions are typically referred to in terms of the "order",
 * which we have denoted via \em L.  In many mathematics texts, the
 * term "degree" is used instead.  Here, both "Legendre order" and
 * "spherical harmonic degree" will refer to the \em l subscript, while
 * the \em m subscript represents the "spherical harmonic order".
 *
 *
 *  @refs
 * - Alain Hebert, <em>Applied Reactor Physics</em>, Presses Internationales
 *   Polytechnique, Montreal, 2009.
 *
 *  @example angle/test/test_SphericalHarmonics.cc
 *
 * Test of SphericalHarmonics.
 *
*/
//---------------------------------------------------------------------------//
class ANGLE_EXPORT SphericalHarmonics
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::size_t size_t;

  //-------------------------------------------------------------------------//
  // SPHERICAL HARMONICS -- a variety of interfaces provided for future work
  //-------------------------------------------------------------------------//

  /**
   *  @brief Calculate \f$ Y^m_l \f$ given \f$\xi\f$ and \f$\varphi\f$.
   *  @param l      Legendre order (or spherical harmonic degree)
   *  @param m      Spherical harmonic order
   *  @param xi     direction cosine w/r to z axis
   *  @param varphi azimutal angle (as defined from x axis)
   *  @return       \f$ Y^m_l(\Omega) \f$
   */
  static double Y_lm(const int l, const int m,
                     const double xi, const double varphi);

  /**
   *  @brief Calculate \f$ Y^m_l \f$ given \f$\mathbf{\Omega}\f$.
   *  @param l   Legendre order (or spherical harmonic degree)
   *  @param m   Spherical harmonic order
   *  @param mu  direction cosine w/r to x axis
   *  @param eta direction cosine w/r to y axis
   *  @param xi  direction cosine w/r to z axis
   *  @return    \f$ Y^m_l(\Omega) \f$
   */
  static double Y_lm(const int l, const int m,
                     const double mu, const double eta, const double xi);

  /**
   *  @brief Calculate \f$ Y^m_l \f$ given \f$\xi\f$.
   *  @param l   Legendre order
   *  @param xi  direction cosine w/r to the polar axis
   *  @return    \f$ Y^0_l(\Omega) = P_l(\xi)\f$
   */
  static double Y_lm(const int l, const double xi);

private:

  /**
   *  @brief Calculate \f$ Y^m_l \f$ given \f$\mathbf{\Omega}\f$.
   *
   *  This private function could be used to interface with
   *  external libraries, leaving the rest of the interfaces
   *  and their implementations unchanged.
   *
   *  @param l   Legendre order (or spherical harmonic degree)
   *  @param m   Spherical harmonic order
   *  @param mu  direction cosine w/r to x axis
   *  @param eta direction cosine w/r to y axis
   *  @param xi  direction cosine w/r to z axis
   *  @return    \f$ Y^m_l(\Omega) \f$
   */
  static double get_Y_lm(const int l, const int m,
                         const double mu, const double eta, const double xi);

};

} // end namespace detran_angle

#endif /* detran_angle_SPHERICAL_HARMONICS_HH_ */

//---------------------------------------------------------------------------//
//              end of SphericalHarmonics.hh
//---------------------------------------------------------------------------//
