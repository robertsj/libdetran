//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TriGauss.hh
 *  @brief TriGauss class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_TRIGAUSS_HH_
#define detran_angle_TRIGAUSS_HH_

#include "angle/BaseQuadrature.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran_angle
{

class ANGLE_EXPORT TriGauss: public BaseQuadrature
{

public:

  static std::string name() {return "tg";}

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    Insist(detran_utilities::soft_equiv(a, -1.0) &&
           detran_utilities::soft_equiv(b,  1.0),
           "Trigauss is defined only for polar angles, [-1, 1].");
    Require(d_m % 2 == 0);
    Require(d_m / 2 <= 6);

    int m = d_m / 2;

    vec_dbl x(m, 0.0);
    vec_dbl w(m, 0.0);

    if (m == 1)
    {
      // not optimal
      x[0] = 0.7071067811865475244008445;
      w[0] = 1.0;
    }
    else if (m == 2)
    {
      // based on [0, pi]
      x[0] = 0.970028106554334;
      x[1] = 0.480049092416837;
      w[0] = 0.144806658914468;
      w[1] = 0.855193341085532;
    }
    else if (m == 3)
    {
      x[0] = 0.983720409069126;
      x[1] = 0.707106781186548;
      x[2] = 0.179705750550371;
      w[0] = 0.079338724686562;
      w[1] = 0.486356182323044;
      w[2] = 0.434305092990394;
    }
    else if (m == 4)
    {
      x[0] = 9.937172857099943e-01;
      x[1] = 8.665969087139184e-01;
      x[2] = 4.990088153604910e-01;
      x[3] = 1.119194178021024e-01;
      w[0] = 3.123105861536087e-02;
      w[1] = 2.526724420559626e-01;
      w[2] = 4.388001783990736e-01;
      w[3] = 2.772963209296030e-01;
    }
    else if (m == 5)
    {
      x[0] = 0.997101416593042;
      x[1] = 0.932984999374027;
      x[2] = 0.707106781186548;
      x[3] = 0.359915255224124;
      x[4] = 0.076083934100109;
      w[0] = 0.014556433410266;
      w[1] = 0.134916240529173;
      w[2] = 0.310026530523421;
      w[3] = 0.349734629912457;
      w[4] = 0.190766165624683;
    }
    else if (m == 6)
    {
      x[0] = 0.998487324664592;
      x[1] = 0.963239238930980;
      x[2] = 0.824493558075628;
      x[3] = 0.565871339344723;
      x[4] = 0.268645060597933;
      x[5] = 0.054982383398188;
      //
      w[0] = 0.007642068937484;
      w[1] = 0.076559758219347;
      w[2] = 0.204518429695601;
      w[3] = 0.297990225104937;
      w[4] = 0.274508539542121;
      w[5] = 0.138780978500511;
    }
    else
    {
      THROW("ERROR!!!");
    }

    for (int i = 0; i < m; ++i)
    {
      d_x[      i] = -x[i];
      d_x[d_m-i-1] =  x[i];
      d_w[      i] =  w[i];
      d_w[d_m-i-1] =  w[i];
    }
  }

};

} // end namespace detran_angle

#endif /* detran_angle_TRIGAUSS_HH_ */

//----------------------------------------------------------------------------//
//              end of file TriGauss.hh
//----------------------------------------------------------------------------//
