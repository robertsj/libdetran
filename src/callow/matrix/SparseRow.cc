//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SparseRow.cc
 *  @brief SparseRow member definitions
 */
//----------------------------------------------------------------------------//

#include "SparseRow.hh"

#include <iostream>
using std::cout;
using std::endl;

namespace callow
{

//----------------------------------------------------------------------------//
SparseRow::SparseRow(const int n)
  : d_n(n)
{
  Require(n >= 0);
  d_elements.reserve(d_n);
}

//----------------------------------------------------------------------------//
SparseRow::SparseRow(const SparseRow& row)
  : d_n(row.size())
  , d_elements(row.elements())
{
  /* ... */
}

//----------------------------------------------------------------------------//
SparseRow::SparseRow(const Matrix &A, const int i)
  : d_n(A.number_rows())
  , d_elements(A.end(i)-A.start(i))
{
  Require(i >= 0);
  Require(i < A.number_rows());

  // copy nonzero elements of A[i, :] into row
  for (int j = 0; j < number_nonzeros(); ++j)
  {
    int k = j + A.start(i);
    d_elements[j].second = A.value(k);
    d_elements[j].first = A.column(k);
  }
}

//----------------------------------------------------------------------------//
SparseRow::SparseRow(const double *a, const int n)
  : d_n(n)
{
  Require(d_n >= 0);

  // copy nonzero elements of A[i, :] into row
  for (int j = 0; j < d_n; ++j)
  {
    if (abs(a[j]) > 0)
    {
      d_elements.push_back(element_t(j, a[j]));
    }
  }
}

//----------------------------------------------------------------------------//
bool SparseRow::verify() const
{
  for (int i = 1; i < d_elements.size(); ++i)
  {
    if (d_elements[i].first <= d_elements[i-1].first)
      return false;
  }
  return true;
}

//----------------------------------------------------------------------------//
bool SparseRow::insert(int i, double v, const int type)
{
  Require(i < d_n);
  if (type == INSERT and v == 0.0)
    return true;
  // find index that contains i
  iterator element = begin();
  for (; element < end(); ++element)
  {
    if (i < element->first)
    {
      break;
    }
    if (i == element->first)
    {
      if (type == INSERT)
      {
        element->second = v;
      }
      else
      {
        double new_v = v + element->second;
        if (new_v == 0)
        {
          d_elements.erase(element);
          return false;
        }
        element->second = new_v;
      }
      return true;
    }
  }
  // as long as indices are ordered, indices[p-1] < i < indices[p]
  d_elements.insert(element, element_t(i, v));
  return true;
}

//----------------------------------------------------------------------------//
bool SparseRow::insert(int a, int b, const SparseRow &r, double v, const int t)
{
  Require(a < b);
  Require(b <= d_n);
  Require(r.size() == d_n);
  Require(r.verify());
  crange_t r_range = r.range(a, b);
  bool rval = true;
  for (auto element = r_range.first; element < r_range.second; ++element)
  {
    rval = insert(element->first, element->second*v, t) || rval;
  }
  return true;
}

//----------------------------------------------------------------------------//
bool SparseRow::update(const double threshold)
{
  return true;
}

//----------------------------------------------------------------------------//
void SparseRow::display(const std::string s) const
{
  printf(" Sparse Row %s \n", s.c_str());
  printf(" ---------------------------\n");
  printf("             size = %5i \n",   d_n);
  printf("  number nonzeros = %5i \n",   number_nonzeros());
  printf("\n");
  for (int p = 0; p < number_nonzeros(); ++p)
  {
      int j = d_elements[p].first;
      double v = d_elements[p].second;
      printf(" %3i (%13.6e)\n", j, v);
  }
  printf("\n");
}

//----------------------------------------------------------------------------//
void SparseRow::clear()
{
  d_elements.clear();
}

//----------------------------------------------------------------------------//
const double& SparseRow::operator[](const int i) const
{
  Require(i >= 0);
  Require(i < d_n);
  int p = find(i);
  Require(i >= 0);
  return d_elements[p].second;
}

//----------------------------------------------------------------------------//
double& SparseRow::operator[](const int i)
{
  Require(i >= 0);
  Require(i < d_n);
  int p = find(i);
  Require(i >= 0);
  return d_elements[p].second;
}

/**
//----------------------------------------------------------------------------//
SparseRow SparseRow::operator()(const int low, const int high)
{
  auto _begin = begin();
  for (; _begin < end(); ++_begin)
  {
    if (_begin->first >= low)  // insertions are before, so we want next one
      break;
  }
  auto _end = _begin+1;
  for (; _end < end(); ++_end)
  {
    if (_end->first >= high) // next one or end
      break;
  }
  return range_t(_begin, _end);
}
*/

//----------------------------------------------------------------------------//
int SparseRow::find(const int i) const
{
  Require(i >= 0);
  Require(i < d_n);
  for (int p = 0; p < number_nonzeros(); ++p)
  {
    if (d_elements[p].first == i)
      return p;
  }
  return -1;
}


//----------------------------------------------------------------------------//
void SparseRow::delete_value(const int p)
{
  Require(p < number_nonzeros());
  d_elements.erase(d_elements.begin()+p);
}

//----------------------------------------------------------------------------//
SparseRow::crange_t SparseRow::range(const int a, const int b) const
{
  auto _begin = cbegin();
  for (; _begin < cend(); ++_begin)
  {
    if (_begin->first >= a)  // insertions are before, so we want next one
      break;
  }
  auto _end = _begin+1;
  for (; _end < cend(); ++_end)
  {
    if (_end->first >= b) // next one or end
      break;
  }
  return crange_t(_begin, _end);
}

//----------------------------------------------------------------------------//
SparseRow::range_t SparseRow::range(const int a, const int b)
{
  auto _begin = begin();
  for (; _begin < cend(); ++_begin)
  {
    if (_begin->first >= a)  // insertions are before, so we want next one
      break;
  }
  auto _end = _begin+1;
  for (; _end < cend(); ++_end)
  {
    if (_end->first > b) // next one or end
      break;
  }
  return range_t(_begin, _end);
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file SparseRow.cc
//----------------------------------------------------------------------------//
