//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCILUT.cc
 *  @brief PCILUT member definitions
 */
//----------------------------------------------------------------------------//

#include "PCILUT.hh"
#include "SparseRow.hh"

#include <cmath>
#include <algorithm>
#include <iostream>
using std::cout;
using std::endl;

#define dbprint(s, i, k, row) void; //printf("%s %i %i ",s, i, k); row.display("");

namespace callow
{


  inline bool drop_t(SparseRow  &r, int idx, double t)
  {
    if (std::abs(r[idx]) < t)
    {
      r.delete_value(idx);
      return false;
    }
    return true;
  }

  bool drop_p(SparseRow  &r, int d, int p, double t)
  {
    // first drop small values
    if (t > 0.0)
    {
      for (int idx = 0; idx < r.number_nonzeros(); ++idx)
      {
        drop_t(r, idx, t);
      }
    }
    // then keep the p largest on either side of d
    int cL = 0;
    int idx_d = 0;
    for (int idx = 0; idx < r.number_nonzeros(); ++idx)
    {
      int j = r.elements()[idx].first;
      if (j == d)
        idx_d = idx;
      else if (j < d)
        cL++;
    }
    int cU = r.number_nonzeros() - 1 - cL;
    Ensure(cU >= 0);

    // copy r and sort the first and second portions
    auto q = r;
    //r.display("r");
    auto rq = q.range(0, d);
    std::sort(rq.first, rq.second, SparseRow::compare_v);
    rq = q.range(d+1, r.size());
    std::sort(rq.first, rq.second, SparseRow::compare_v);
    // iterators are invalid now.  however, we know we need to keep
    // the largest p (or drop the smallest)
    int drop = cL - p;
    if (drop > 0)
    {
    	q.elements().erase(q.begin(), q.begin()+drop);
    	cL = p;
    }
    drop = cU - p;
    if (drop>0)
    {
    	q.elements().erase(q.begin()+cL+1, q.begin()+cL+1+drop);
    	cU = p;
    }
    // sort q by column and copy to r
    std::sort(q.begin(), q.end(), SparseRow::compare_c);
    r = q;
    return true;
  }

//----------------------------------------------------------------------------//
PCILUT::PCILUT(SP_matrix A, const size_t p, const double t)
  : Base("PCILUT")
  , d_p(p)
  , d_t(t)
{
  Require(A);
  Require(A->number_rows() == A->number_columns());
  Insist(dynamic_cast<Matrix*>(A.bp()),
    "Need an explicit matrix for use with PCILUT");


  //printf("***** PCILUT(%i, %6.2f)\n", d_p, d_t);

  SP_matrixfull B(A);

  d_size = A->number_columns();

  // get the csr structure
  int n = B->number_rows();

  std::vector<SparseRow> rows;
  bool keep = false;
  for (int i = 0; i < n; ++i)
  {
    SparseRow row(*B, i);
    for (int k = 0; k < i; ++k)
    {
      int p = row.find(k);
      if (p < 0)
      {
        continue;
      }
      SparseRow &row_U = rows[k];
      double diag_val = row_U.elements()[row_U.find(k)].second;
      Ensure(diag_val != 0);
      row.elements()[p].second /= diag_val;
      //keep = false;
      //if (d_t > 0)
      keep = drop_t(row, p, d_t);
      if (!keep)
        continue;
      dbprint("-->", i, k, row);
      row.insert(k+1, n, row_U, -row.elements()[p].second, SparseRow::ADD);
      dbprint("<--", i, k, row);

    } // end k
    keep = drop_p(row, i, d_p, d_t);
    //row.display("final");
    rows.push_back(row);
  }

  // copy A
  d_P = new Matrix(n, n, 2*d_p+1);
  for (int i = 0; i < n; ++i)
  {
    for (auto it = rows[i].begin(); it < rows[i].end(); ++it)
    {
      d_P->insert(i, it->first, it->second);
	  }
  }
  d_P->assemble();
  printf("PCILUT...construction complete.\n");
  // d_P->print_matlab("pcilut.out");
}

//----------------------------------------------------------------------------//
PCILUT::SP_preconditioner
PCILUT::Create(SP_matrix A, const size_t p, const double t)
{
  SP_preconditioner pc(new PCILUT(A, p, t));
  return pc;
}

//----------------------------------------------------------------------------//
void PCILUT::display(const std::string &name)
{
  cout << " PCILU " << name << endl;
  d_P->display();
}

//----------------------------------------------------------------------------//
void PCILUT::apply(Vector &b, Vector &x)
{
	// solve LUx = x --> x = inv(U)*inv(L)*x

	// forward substitution
	//   for i = 0:m-1
	//     x[i] = 1/L[i,i] * ( b[i] - sum(k=0:i-1, L[i,k]*y[k]) )
	// but note that in our ILU(0) scheme, L is *unit* lower triangle,
	// meaning L has ones on the diagonal (whereas U does not)
	Vector y(b.size(), 0.0);

	for (int i = 0; i < d_P->number_rows(); ++i)
	{
		// start index
		int s = d_P->start(i);
		// diagonal index
		int d = d_P->diagonal(i);
		// invert row
		y[i] = b[i];
		for (int p = s; p < d; ++p)
		{
			// column index
			int c = d_P->column(p);
			y[i] -= d_P->values()[p] * y[c];
		}
	}

	// backward substitution
	//   for i = m-1:0
	//     y[i] = 1/U[i,i] * ( b[i] - sum(k=i+1:m-1, U[i,k]*y[k]) )
	for (int i = d_P->number_rows() - 1; i >= 0; --i)
	{
		// diagonal index
		int d = d_P->diagonal(i);
		// end index
		int e = d_P->end(i);
		// invert row
		x[i] = y[i];
		for (int p = d+1; p < e; ++p)
		{
			// column index
			int c = d_P->column(p);
			x[i] -= d_P->values()[p] * x[c];
		}
		x[i] /= d_P->values()[d];
	}
}


} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file PCILU0.cc
//----------------------------------------------------------------------------//
