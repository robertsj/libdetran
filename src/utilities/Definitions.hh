//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Definitions.hh
 *  @author Jeremy Roberts
 *  @brief  Useful type definitions, etc.
 */
//---------------------------------------------------------------------------//

#ifndef DETRAN_UTILS_DEFINITIONS_HH_
#define DETRAN_UTILS_DEFINITIONS_HH_

#include "utilities/utilities_export.hh"
#include <vector>
#include <map>

namespace detran_utilities
{

typedef unsigned int size_t;
typedef unsigned int u_int;

// Templated standard vec

template <typename T>
struct vec1_T {typedef std::vector<T> value_type;};
template <typename T>
struct vec2_T {typedef std::vector<typename vec1_T<T>::value_type> value_type;};
template <typename T>
struct vec3_T {typedef std::vector<typename vec2_T<T>::value_type> value_type;};
template <typename T>
struct vec4_T {typedef std::vector<typename vec3_T<T>::value_type> value_type;};

// Instantiations

typedef std::vector<int> vec_int;
typedef std::vector<std::vector<int> > vec2_int;
typedef std::vector<std::vector<std::vector<int> > > vec3_int;

typedef std::vector<double> vec_dbl;
typedef std::vector<std::vector<double> > vec2_dbl;
typedef std::vector<std::vector<std::vector<double> > > vec3_dbl;
typedef std::vector<std::vector<std::vector<std::vector<double> > > > vec4_dbl;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > vec5_dbl;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > > vec6_dbl;

typedef std::vector<size_t> vec_size_t;
typedef std::vector<std::vector<size_t> > vec2_size_t;
typedef std::vector<std::vector<std::vector<size_t> > > vec3_size_t;
typedef std::vector<std::vector<std::vector<std::vector<size_t> > > > vec4_size_t;

typedef std::vector<bool> vec_bool;

// DLL exports
UTILITIES_TEMPLATE_EXPORT(std::vector<int>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec_int>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec2_int>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec3_int>)
UTILITIES_TEMPLATE_EXPORT(std::vector<size_t>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec_size_t>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec2_size_t>)
UTILITIES_TEMPLATE_EXPORT(std::vector<double>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec_dbl>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec2_dbl>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec3_dbl>)
UTILITIES_TEMPLATE_EXPORT(std::vector<vec4_dbl>)
UTILITIES_TEMPLATE_EXPORT(std::basic_string<char, std::char_traits<char>, std::allocator<char> >)
UTILITIES_TEMPLATE_EXPORT(std::map<std::string, int>)
UTILITIES_TEMPLATE_EXPORT(std::map<std::string, double>)
UTILITIES_TEMPLATE_EXPORT(std::map<std::string, vec_int>)
UTILITIES_TEMPLATE_EXPORT(std::map<std::string, vec_dbl>)
UTILITIES_TEMPLATE_EXPORT(std::map<std::string, std::string>)
UTILITIES_TEMPLATE_EXPORT(std::vector<std::_Vbase, std::allocator<bool> >)
UTILITIES_TEMPLATE_EXPORT(std::vector<bool>)

//disable warnings on 255 char debug symbols
#pragma warning (disable : 4244)

} // end namespace detran_utilities

#endif // DETRAN_UTILS_DEFINITIONS_HH_
