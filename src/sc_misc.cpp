// ===============================================================
//
// SCArray R package
// Copyright (C) 2023    Xiuwen Zheng
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "sc_main.h"


extern "C"
{
// ====  rowCollapse & colCollapse  ====

static int idx_n = 0, idx_i = 0, idx_st = 0;
static int *idx_p = NULL;

SEXP c_Collapse_init(SEXP idx)
{
	idx_n = Rf_length(idx);
	idx_i = 0; idx_st = 0;
	idx_p = (idx_n > 0) ? INTEGER(idx) : &NA_INTEGER;
	return R_NilValue;
}


/*
SEXP c_rowCollapse(SEXP mat, SEXP val, SEXP idx)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	double *pv = REAL(val);
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i_j {
			double v = p[j], &d = pv[j];
			if (!ISNAN(v))
			{
				if (!ISNAN(d) && (v < d)) d = v;
			} else {
				if (!na_rm) d = NA_REAL;
			}			
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i_j {
			int v = p[j]; double &d = pv[j];
			if (v != NA_INTEGER)
			{
				if (!ISNAN(d) && (v < d)) d = v;
			} else {
				if (!na_rm) d = NA_REAL;
			}			
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(nrow)));
		memset(pn, 0, sizeof(int)*nrow);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double v = p[i], &d = pv[r];
				if (!ISNAN(v))
				{
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double &d = pv[r];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < nrow; i++)
			if (pn[i]<ncol && !ISNAN(pv[i]) && 0<pv[i]) pv[i] = 0;
		UNPROTECT(1);
	}
	// output
	return val;
}
*/

SEXP c_colCollapse(SEXP mat, SEXP val)
{
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat), *pv = REAL(val);
		FOR_LOOP_i {
			int j = idx_p[idx_i++];
			if (idx_i >= idx_n) idx_i = 0;
			pv[idx_st+i] = (j != NA_INTEGER) ? p[j-1] : NA_REAL;
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		int *p = INTEGER(mat), *pv = INTEGER(val);
		FOR_LOOP_i {
			int j = idx_p[idx_i++];
			if (idx_i >= idx_n) idx_i = 0;
			pv[idx_st+i] = (j != NA_INTEGER) ? p[j-1] : j;
		}
	}
	idx_st += ncol;
	// output
	return val;
}

}
