//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ColorMap.cc
 *  @brief ColorMap member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ColorMap.hh"
#include  <cstdio>

namespace detran_ioutils
{

//----------------------------------------------------------------------------//
ColorMap::vec_rgb
ColorMap::color(const size_t map, const vec_dbl &values)
{
  Require(map < END_COLOR_MAPS);

  using namespace detran_utilities;

  if (map == DEFAULT or map == COOLWARM)
  {
    return colors_function(map, values);
  }
  else if (map == RANDOM)
  {
    return colors_random(values);
  }
  else
  {
    vec_dbl scaled = values;
    double min_val = vec_min(scaled);
    double max_val = vec_max(scaled);
    vec_plus_a(scaled, -min_val);
    if (max_val - min_val > 0)
      vec_scale(scaled, 1.0 / (max_val - min_val));
    return colors_interp(map, scaled);
  }
}

//----------------------------------------------------------------------------//
ColorMap::vec_rgb
ColorMap::colors_function(const size_t map, const vec_dbl &values)
{
  ColorMap::vec_rgb c(values.size());

  double min = detran_utilities::vec_min(values);
  double max = detran_utilities::vec_max(values);

  for (size_t i = 0; i < values.size(); ++i)
  {
    unsigned char r = 255;
    unsigned char g = 255;
    unsigned char b = 255;
    if (values[i] >= 0.0)
    {
      double v = (values[i]-min)/(max-min);
      if (map == DEFAULT)
      {
        c[i].r = 255. * (2.0*v - v*v);
        c[i].g = 255. * (4.0*v - 4.0*v*v);
        c[i].b = 255. * (1.0-v*v);
      }
      else if (map == COOLWARM)
      {
        c[i].r = 255. * (0.5764 - 0.3312*cos(3.788*v) + 0.2042*sin(3.788*v));
        c[i].g = 255. * (0.3989 - 0.1124*cos(3.914*v) + 0.4392*sin(3.914*v));
        c[i].b = 255. * (0.5932 + 0.1619*cos(4.024*v) + 0.3747*sin(4.024*v));
      }
    }
  }
  return c;
}

//----------------------------------------------------------------------------//
ColorMap::vec_rgb
ColorMap::colors_random(const vec_dbl &values)
{
  using std::find;

  typedef std::vector<gray_to_color> vec_g2c;

  // unique values
  vec_dbl unique = detran_utilities::vec_unique(values);

  // value to color map
  vec_g2c unique_c(unique.size());
  detran_utilities::Random R(1234);
  for (size_t i = 0; i < unique_c.size(); i++)
  {
    unique_c[i].value   = unique[i];
    unique_c[i].color.r = 255. * R.rnd();
    unique_c[i].color.g = 255. * R.rnd();
    unique_c[i].color.b = 255. * R.rnd();
  }

  // colors to return
  vec_rgb c(values.size());
  for (size_t i = 0; i < values.size(); ++i)
  {
    vec_g2c::iterator it = find(unique_c.begin(), unique_c.end(), values[i]);
    Assert(it != unique_c.end());
    c[i] = it->color;
  }

  return c;
}

//----------------------------------------------------------------------------//
ColorMap::vec_rgb
ColorMap::colors_interp(const size_t map, const vec_dbl &values)
{
  using std::sqrt;
  using std::pow;

  // color assignments
  int jet[9] =
  { darkblue, blue, azure, cyan, lightgreen, yellow, orange, red, darkred };
  int hot[4] =
  { black, red, yellow, white };
  int cool[5] =
  { cyan, fuchsia };

  int *cmap;
  int n;
  if (map == JET)
  {
    cmap = jet;
    n = 9;
  }
  else if (map == HOT)
  {
    cmap = hot;
    n = 4;
  }
  else if (map == COOL)
  {
    cmap = cool;
    n = 2;
  }

  // color node edges in 1-D
  vec_dbl d(n, 0.0);
  for (int i = 1; i < n; ++i)
  {
    rgb_t c0 = hex_to_rgb(cmap[i - 1]);
    rgb_t c1 = hex_to_rgb(cmap[i    ]);
    d[i] = d[i-1] +
           sqrt(pow(c1.r-c0.r, 2) + pow(c1.b-c0.b, 2) + pow(c1.g-c0.g, 2));
  }
  double total_d = d[n-1];
  detran_utilities::vec_scale(d, 1.0 / total_d);

  // loop through all values and find the appropriate bin and interpolate
  ColorMap::vec_rgb c(values.size());
  for (int i = 0; i < values.size(); ++i)
  {
    int j = 0;
    for (; j < d.size(); ++j) if (d[j] >= values[i]) break;
    if (j == 0)
    {
      c[i] = hex_to_rgb(cmap[0]);
    }
    else if (j == d.size())
    {
      c[i] = hex_to_rgb(cmap[n-1]);
    }
    else
    {
      rgb_t c0 = hex_to_rgb(cmap[j-1]);
      rgb_t c1 = hex_to_rgb(cmap[j  ]);
      double del = values[i] - d[j-1];
      c[i].r = c0.r + (c1.r - c0.r) * del / (d[j] - d[j-1]);
      c[i].g = c0.g + (c1.g - c0.g) * del / (d[j] - d[j-1]);
      c[i].b = c0.b + (c1.b - c0.b) * del / (d[j] - d[j-1]);
    }
  }
  return c;
}

//----------------------------------------------------------------------------//
ColorMap::rgb_t ColorMap::hex_to_rgb(const int hex)
{
  ColorMap::rgb_t c;
  c.r = (hex >> 16) & 0xFF;
  c.g = (hex >> 8)  & 0xFF;
  c.b = hex & 0xFF;
  return c;
}

} // end namespace detran_ioutils

//----------------------------------------------------------------------------//
//              end of file ColorMap.cc
//----------------------------------------------------------------------------//
