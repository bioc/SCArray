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

#include "sc_main.h"


// SparseArraySeed
SparseMatrix::SparseMatrix(SEXP mat)
{
	static const char *err_type =
		"SparseArraySeed should be a numeric matrix.";
	// slot nzindex
	SEXP ii = GET_SLOT(mat, mkString("nzindex"));
	if (!Rf_isMatrix(ii)) Rf_error(err_type);
	nnzero = INTEGER(GET_DIM(ii))[0];
	nzi_r = INTEGER(ii); nzi_c = nzi_r + nnzero;
	// slot nzdata
	nzdata = GET_SLOT(mat, mkString("nzdata"));
	if (TYPEOF(nzdata)!=REALSXP && TYPEOF(nzdata)!=INTSXP)
		Rf_error(err_type);
	if (Rf_length(nzdata) != nnzero)
		Rf_error(err_type);
}


static void throw_error_type(SEXP x)
{
	SEXP nm = Rf_getAttrib(x, R_ClassSymbol);
	if (Rf_isString(nm))
	{
		const char *s = CHAR(STRING_ELT(nm, 0));
		Rf_error("Unsupported data type 'class: %s'.", s);
	} else {
		Rf_error("Unsupported data type 'typeof: %d'.", (int)TYPEOF(x));
	}
}


bool is_sparse_seed(SEXP mat)
{
	return Rf_inherits(mat, "SparseArraySeed") == TRUE;
}

void get_mat_size(SEXP mat, int &nrow, int &ncol)
{
	static const char *err = "The input should be a numeric matrix.";
	if (Rf_isMatrix(mat))
	{
		int *p = INTEGER(GET_DIM(mat));
		nrow = p[0]; ncol = p[1];
		if (TYPEOF(mat)!=REALSXP && TYPEOF(mat)!=INTSXP)
			throw_error_type(mat);
	} else if (is_sparse_seed(mat))
	{
		SEXP dm = GET_SLOT(mat, mkString("dim"));
		if (Rf_isNull(dm) || Rf_length(dm) != 2) Rf_error(err);
		int *p = INTEGER(dm);
		nrow = p[0]; ncol = p[1];
	} else
		throw_error_type(mat);
}



extern "C"
{

// ====  rowSums  ====

SEXP c_rowSums(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pv = REAL(val);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) pv[j] += p[j];
		} else {
			FOR_LOOP_i_j  pv[j] += p[j];
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) pv[j] += p[j];
		} else {
			FOR_LOOP_i_j
				if (p[j] != NA_INTEGER) pv[j] += p[j]; else { pv[j] = NA_REAL; break; }
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  if (!ISNAN(p[i])) pv[M.nzi_r[i]-1] += p[i];
			} else {
				FOR_LOOP_Mi  pv[M.nzi_r[i]-1] += p[i];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  if (p[i] != NA_INTEGER) pv[M.nzi_r[i]-1] += p[i];
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_r[i]-1] += p[i];
					else pv[M.nzi_r[i]-1] = NA_REAL;
			}
		}
	}
	// output
	return val;
}


// ====  rowProds  ====

SEXP c_rowProds(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	// do
	double *pv = REAL(val);
	if (TYPEOF(mat) == REALSXP)
	{
		get_mat_size(mat, nrow, ncol);
		double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) pv[j] *= p[j];
		} else {
			FOR_LOOP_i_j  pv[j] *= p[j];
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		get_mat_size(mat, nrow, ncol);
		int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) pv[j] *= p[j];
		} else {
			FOR_LOOP_i_j
				if (p[j] != NA_INTEGER) pv[j] *= p[j]; else pv[j] = NA_REAL;
		}
	} else 
		throw_error_type(mat);
	// output
	return val;
}


// ====  rowMeans & colMeans  ====

SEXP c_rowMeans(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pv = REAL(val);
	int *pn = (int*)&pv[nrow];
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) { pv[j] += p[j]; pn[j]++; }
		} else {
			FOR_LOOP_i_j  pv[j] += p[j];
			for (int j=0; j < nrow; j++) pn[j] += ncol;
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) { pv[j] += p[j]; pn[j]++; }
		} else {
			FOR_LOOP_i_j
				if (p[j] != NA_INTEGER) pv[j] += p[j]; else { pv[j] = NA_REAL; break; }
			for (int j=0; j < nrow; j++) pn[j] += ncol;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  if (!ISNAN(p[i]))
					{ int r=M.nzi_r[i]-1; pv[r] += p[i]; pn[r]++; }
			} else {
				FOR_LOOP_Mi  pv[M.nzi_r[i]-1] += p[i];
				for (int j=0; j < nrow; j++) pn[j] += ncol;
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  if (p[i] != NA_INTEGER)
					{ int r=M.nzi_r[i]-1; pv[r] += p[i]; pn[r]++; }
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_r[i]-1] += p[i];
					else pv[M.nzi_r[i]-1] = NA_REAL;
				for (int j=0; j < nrow; j++) pn[j] += ncol;
			}
		}
	}
	// output
	return val;
}

SEXP c_rowMeans_final(SEXP val)
{
	const size_t n = Rf_xlength(val) / 2;
	const double *sv = REAL(val);
	const int *sn = (int*)&sv[n];
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	for (size_t i=0; i < n; i++) p[i] = sv[i] / sn[i];
	return ans;
}

SEXP c_colMeans(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Avg = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Avg);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0; int n=0;
				for (int j=0; j < nrow; j++)
					if (!ISNAN(p[j])) { s += p[j]; n++; }
				pv[i] = s / n;
			}
		} else {
			FOR_LOOP_i {
				double s=0;
				for (int j=0; j < nrow; j++) s += p[j];
				pv[i] = s / nrow;
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0; int n=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) { s += p[j]; n++; }
				pv[i] = s / n;
			}
		} else {
			FOR_LOOP_i {
				double s=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) s += p[j]; else { s = NA_REAL; break; }
				pv[i] = s / nrow;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		memset(pv, 0, sizeof(double)*ncol);
		SEXP Num = PROTECT(NEW_INTEGER(ncol));
		int *pn = INTEGER(Num);
		for (int i=0; i < ncol; i++) pn[i] = nrow;
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  if (!ISNAN(p[i]))
					pv[M.nzi_c[i]-1] += p[i]; else pn[M.nzi_c[i]-1]--;
			} else {
				FOR_LOOP_Mi  pv[M.nzi_c[i]-1] += p[i];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  if (p[i] != NA_INTEGER)
					pv[M.nzi_c[i]-1] += p[i]; else pn[M.nzi_c[i]-1]--;
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_c[i]-1] += p[i];
					else pv[M.nzi_c[i]-1] = NA_REAL;
			}
		}
		for (int i=0; i < ncol; i++) pv[i] /= pn[i];
		UNPROTECT(1);
	}
	// output
	UNPROTECT(1);
	return Avg;
}


// ====  rowVars  ====

SEXP c_rowVars(SEXP mat, SEXP val, SEXP narm, SEXP center)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const bool no_center = (Rf_isNull(center) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pS = REAL(val), *pS2 = pS+nrow;
	int *pN = (int*)&pS[2*nrow];
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat), v;
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i_j  if (!ISNAN(v=p[j]))
					{ pS[j] += v; pS2[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  { v = p[j]; pS[j] += v; pS2[j] += v*v; }
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i_j  if (!ISNAN(p[j]))
					{ v = p[j]-pC[j]; pS[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  { v = p[j]-pC[j]; pS[j] += v*v; }
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		int *p = INTEGER(mat); double v;
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]; pS[j] += v; pS2[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]; pS[j] += v; pS2[j] += v*v; }
					else { pS[j] = pS2[j] = NA_REAL; break; }
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]-pC[j]; pS[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]-pC[j]; pS[j] += v*v; } else { pS[j] = NA_REAL; break; }
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		if (na_rm)
			for (int j=0; j < nrow; j++) pN[j] += ncol;
		if (!no_center)
		{
			const double *pC = REAL(center);
			for (int j=0; j < nrow; j++) pS[j] += ncol*pC[j]*pC[j];
		}
		int r;
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata), v;
			if (no_center)
			{
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (!ISNAN(v=p[i])) { pS[r] += v; pS2[r] += v*v; } else pN[r]--;
					}
				} else {
					FOR_LOOP_Mi
						{ r = M.nzi_r[i]-1; v = p[i]; pS[r] += v; pS2[r] += v*v; }
				}
			} else {
				const double *pC = REAL(center);
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (!ISNAN(v=p[i])) pS[r] += v*v - 2*v*pC[r];
						else { pN[r]--; pS[r] -= pC[r]*pC[r]; }
					}
				} else {
					FOR_LOOP_Mi
						{ r = M.nzi_r[i]-1; v = p[i]; pS[r] += v*v - 2*v*pC[r]; }
				}
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP
			double v;
			if (no_center)
			{
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER) { v = p[i]; pS[r] += v; pS2[r] += v*v; } else pN[r]--;
					}
				} else {
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pS[r] += v; pS2[r] += v*v; }
						else { pS[r] = pS2[r] = NA_REAL; }
					}
				}
			} else {
				const double *pC = REAL(center);
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER) { v = p[i]; pS[r] += v*v - 2*v*pC[r]; }
						else { pN[r]--; pS[r] -= pC[r]*pC[r]; }
					}
				} else {
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER) { v = p[i]; pS[r] += v*v - 2*v*pC[r]; }
						else pS[r] = NA_REAL;
					}
				}
			}
		}
	}
	// output
	if (!na_rm) for (int j=0; j < nrow; j++) pN[j] += ncol;
	return val;
}

SEXP c_rowVars_final(SEXP val, SEXP center)
{
	const bool no_center = Rf_isNull(center) == TRUE;
	const size_t n = Rf_xlength(val) / 3;
	const double *pS = REAL(val), *pS2 = pS+n;
	const int *pN = (const int*)&pS[2*n];
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	if (no_center)
	{
		for (size_t i=0; i < n; i++)
			p[i] = (pN[i]>1) ? (pS2[i] - pS[i]*pS[i]/pN[i]) / (pN[i]-1) : NA_REAL;
	} else {
		for (size_t i=0; i < n; i++)
			p[i] = (pN[i]>1) ? pS[i] / (pN[i]-1) : NA_REAL;
	}
	return ans;
}

SEXP c_colVars(SEXP mat, SEXP narm, SEXP center)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const bool no_center = (Rf_isNull(center) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Var = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Var);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, s2=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (!ISNAN(v=p[j])) { s += v; s2 += v*v; n++; }
					pv[i] = (n>1) ? (s2 - s*s/n) / (n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, s2=0, v;
					for (int j=0; j < nrow; j++) { v = p[j]; s += v; s2 += v*v; }
					pv[i] = (nrow>1) ? (s2 - s*s/nrow) / (nrow-1) : NA_REAL;
				}
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (!ISNAN(p[j])) { v = p[j]-pC[i]; s += v*v; n++; }
					pv[i] = (n>1) ? s/(n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, v;
					for (int j=0; j < nrow; j++) { v = p[j]-pC[i]; s += v*v; }
					pv[i] = (nrow>1) ? s/(nrow-1) : NA_REAL;
				}
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, s2=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]; s += v; s2 += v*v; n++; }
					pv[i] = (n>1) ? (s2 - s*s/n) / (n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, s2=0, v;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]; s += v; s2 += v*v; }
						else { s = s2 = NA_REAL; break; }
					pv[i] = (nrow>1) ? (s2 - s*s/nrow) / (nrow-1) : NA_REAL;
				}
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]-pC[i]; s += v*v; n++; }
					pv[i] = (n>1) ? s/(n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, v;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]-pC[i]; s += v*v; }
						else { s = NA_REAL; break; }
					pv[i] = (nrow>1) ? s/(nrow-1) : NA_REAL;
				}
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		int *pN=NULL;
		if (na_rm)
		{
			pN = INTEGER(PROTECT(NEW_INTEGER(ncol)));
			for (int i=0; i < ncol; i++) pN[i] = nrow;
		}
		// do sparse
		if (no_center)
		{
			memset(pv, 0, sizeof(double)*ncol);
			double *p2 = REAL(PROTECT(NEW_NUMERIC(ncol)));
			memset(p2, 0, sizeof(double)*ncol);
			if (TYPEOF(M.nzdata) == REALSXP)
			{
				double *p = REAL(M.nzdata), v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (!ISNAN(v=p[i])) { pv[c] += v; p2[c] += v*v; } else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						v = p[i]; pv[c] += v; p2[c] += v*v;
					}
				}
			} else {
				int *p = INTEGER(M.nzdata);  // INTSXP
				double v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v; p2[c] += v*v; } else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v; p2[c] += v*v; } else pv[c] = p2[c] = NA_REAL;
					}
				}
			}
			// final
			for (int i=0; i < ncol; i++)
			{
				int n = na_rm ? pN[i] : nrow;
				pv[i] = (n>1) ? (p2[i] - pv[i]*pv[i]/n) / (n-1) : NA_REAL;
			}
		} else {
			const double *pC = REAL(center);
			for (int i=0; i < ncol; i++) pv[i] = nrow*pC[i]*pC[i];
			if (TYPEOF(M.nzdata) == REALSXP)
			{
				double *p = REAL(M.nzdata), v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (!ISNAN(v=p[i])) pv[c] += v*v - 2*v*pC[c]; else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						v = p[i]; pv[c] += v*v - 2*v*pC[c];
					}
				}
			} else {
				int *p = INTEGER(M.nzdata);  // INTSXP
				double v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v*v - 2*v*pC[c]; } else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v*v - 2*v*pC[c]; } else pv[c] = NA_REAL;
					}
				}
			}
			// final
			for (int i=0; i < ncol; i++)
			{
				int n = na_rm ? pN[i] : nrow;
				pv[i] = (n>1) ? pv[i]/(n-1) : NA_REAL;
			}
		}
		// release
		UNPROTECT(no_center + na_rm);
	}
	// output
	UNPROTECT(1);
	return Var;
}


}
