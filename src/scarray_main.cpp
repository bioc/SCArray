// ===============================================================
//
// SCArray R package
// Copyright (C) 2021-2023    Xiuwen Zheng
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

#include <R.h>
#include <Rdefines.h>


extern "C"
{

inline static void get_nrow_ncol(SEXP mat, int &nrow, int &ncol)
{
	SEXP dm = GET_DIM(mat);
	if (Rf_isNull(dm) || Rf_length(dm) != 2)
		Rf_error("It should be a numeric matrix.");
	int *pdm = INTEGER(dm);
	nrow = pdm[0];
	ncol = pdm[1];
}

// ====  rowSums  ====

SEXP c_rowSums_update(SEXP mat, SEXP val, SEXP narm)
{
	// get nrow and ncol
	int nrow, ncol;
	get_nrow_ncol(mat, nrow, ncol);
	// do
	const double *p = REAL(mat);
	double *pv = REAL(val), v;
	if (Rf_asLogical(narm) == TRUE)
	{
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
				if (!ISNAN(v=p[j])) pv[j] += v;
		}
	} else {
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
				pv[j] += p[j];
		}
	}
	// output
	return val;
}


// ====  rowProds  ====

SEXP c_rowProds_update(SEXP mat, SEXP val, SEXP narm)
{
	// get nrow and ncol
	int nrow, ncol;
	get_nrow_ncol(mat, nrow, ncol);
	// do
	const double *p = REAL(mat);
	double *pv = REAL(val), v;
	if (Rf_asLogical(narm) == TRUE)
	{
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
				if (!ISNAN(v=p[j])) pv[j] *= v;
		}
	} else {
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
				pv[j] *= p[j];
		}
	}
	// output
	return val;
}


// ====  rowMeans  ====

SEXP c_rowMeans_update(SEXP mat, SEXP val, SEXP narm)
{
	// get nrow and ncol
	int nrow, ncol;
	get_nrow_ncol(mat, nrow, ncol);
	// do
	const double *p = REAL(mat);
	double *pv = REAL(val), *pn = pv+nrow, v;
	if (Rf_asLogical(narm) == TRUE)
	{
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
			{
				if (!ISNAN(v=p[j]))
				{
					pv[j] += v; pn[j] ++;
				}
			}
		}
	} else {
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
				pv[j] += p[j];
		}
		for (int j=0; j < nrow; j++) pn[j] += ncol;
	}
	// output
	return val;
}

SEXP c_rowMeans_final(SEXP val)
{
	const size_t n = Rf_xlength(val) / 2;
	const double *sv = REAL(val), *sn = sv + n;
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	for (size_t i=0; i < n; i++) p[i] = sv[i] / sn[i];
	return ans;
}


// ====  rowVars  ====

SEXP c_rowVars_update(SEXP mat, SEXP val, SEXP narm, SEXP center)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	// get nrow and ncol
	int nrow, ncol;
	get_nrow_ncol(mat, nrow, ncol);
	// do
	double *pS=REAL(val), *pS2=pS+nrow, *pN=pS+2*nrow, v;
	const double *p = REAL(mat);
	if (Rf_isNull(center))
	{
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
			{
				v = p[j];
				if (!na_rm || !ISNAN(v))
				{
					pS[j] += v; pS2[j] += v*v; pN[j] ++;
				}
			}
		}
	} else {
		const double *pC = REAL(center);
		for (int i=0; i < ncol; i++, p+=nrow)
		{
			for (int j=0; j < nrow; j++)
			{
				if (!na_rm || !ISNAN(p[j]))
				{
					v = p[j] - pC[j]; pS[j] += v*v; pN[j] ++;
				}
			}
		}
	}
	// output
	return val;
}

SEXP c_rowVars_final(SEXP val, SEXP center)
{
	const bool has_no_center = Rf_isNull(center);
	const size_t n = Rf_xlength(val) / 3;
	const double *pS=REAL(val), *pS2=pS+n, *pN=pS+2*n;
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	if (has_no_center)
	{
		for (size_t i=0; i < n; i++)
			p[i] = (pN[i]>1) ? (pS2[i] - pS[i]*pS[i]/pN[i]) / (pN[i]-1) : NA_REAL;
	} else {
		for (size_t i=0; i < n; i++)
			p[i] = (pN[i]>1) ? pS[i] / (pN[i]-1) : NA_REAL;
	}
	return ans;
}


}
