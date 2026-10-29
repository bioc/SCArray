// ===============================================================
//
// SCArray R package
// Copyright (C) 2024-2026    Xiuwen Zheng
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

// Registration of native routines, so that .Call() checks the number of
// arguments and symbols are not searched dynamically.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

// The C++ kernels are declared 'extern "C"', so they have C linkage here.

extern SEXP c_add(SEXP, SEXP);
extern SEXP c_get_split(SEXP, SEXP, SEXP);
extern SEXP c_sparse_blocksize(SEXP, SEXP, SEXP, SEXP);

extern SEXP c_rowSums(SEXP, SEXP, SEXP);
extern SEXP c_colSums(SEXP, SEXP);
extern SEXP c_rowProds(SEXP, SEXP, SEXP);
extern SEXP c_rowProds_final(SEXP);
extern SEXP c_colProds(SEXP, SEXP);

extern SEXP c_rowMeans(SEXP, SEXP, SEXP);
extern SEXP c_rowMeans_final(SEXP);
extern SEXP c_colMeans(SEXP, SEXP);
extern SEXP c_rowWMeans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_rowWMeans_final(SEXP);
extern SEXP c_colWMeans(SEXP, SEXP, SEXP);

extern SEXP c_rowVars(SEXP, SEXP, SEXP, SEXP);
extern SEXP c_rowVars_final(SEXP, SEXP);
extern SEXP c_colVars(SEXP, SEXP, SEXP);
extern SEXP c_rowWVars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_rowWVars_final(SEXP);
extern SEXP c_colWVars(SEXP, SEXP, SEXP);

extern SEXP c_rowMeanVar_final(SEXP);
extern SEXP c_colMeanVar(SEXP, SEXP);

extern SEXP c_rowMins(SEXP, SEXP, SEXP);
extern SEXP c_colMins(SEXP, SEXP);
extern SEXP c_rowMaxs(SEXP, SEXP, SEXP);
extern SEXP c_colMaxs(SEXP, SEXP);
extern SEXP c_rowRanges(SEXP, SEXP, SEXP);
extern SEXP c_colRanges(SEXP, SEXP);

extern SEXP c_row_sum_grp(SEXP, SEXP, SEXP, SEXP);
extern SEXP c_col_sum_grp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_row_nnzero(SEXP, SEXP, SEXP);
extern SEXP c_col_nnzero(SEXP, SEXP);
extern SEXP c_row_anyNA(SEXP, SEXP);
extern SEXP c_col_anyNA(SEXP);

extern SEXP c_rowCollapse_init(SEXP, SEXP);
extern SEXP c_rowCollapse(SEXP, SEXP, SEXP, SEXP);
extern SEXP c_colCollapse(SEXP, SEXP, SEXP, SEXP);


#define CALL(name, nargs)    {#name, (DL_FUNC)&name, nargs}

static const R_CallMethodDef CallEntries[] =
{
	CALL(c_add, 2),
	CALL(c_get_split, 3),
	CALL(c_sparse_blocksize, 4),

	CALL(c_rowSums, 3),   CALL(c_colSums, 2),
	CALL(c_rowProds, 3),  CALL(c_rowProds_final, 1),  CALL(c_colProds, 2),

	CALL(c_rowMeans, 3),  CALL(c_rowMeans_final, 1),  CALL(c_colMeans, 2),
	CALL(c_rowWMeans, 6), CALL(c_rowWMeans_final, 1), CALL(c_colWMeans, 3),

	CALL(c_rowVars, 4),  CALL(c_rowVars_final, 2),  CALL(c_colVars, 3),
	CALL(c_rowWVars, 6), CALL(c_rowWVars_final, 1), CALL(c_colWVars, 3),

	CALL(c_rowMeanVar_final, 1),  CALL(c_colMeanVar, 2),

	CALL(c_rowMins, 3),    CALL(c_colMins, 2),
	CALL(c_rowMaxs, 3),    CALL(c_colMaxs, 2),
	CALL(c_rowRanges, 3),  CALL(c_colRanges, 2),

	CALL(c_row_sum_grp, 4),  CALL(c_col_sum_grp, 6),
	CALL(c_row_nnzero, 3),   CALL(c_col_nnzero, 2),
	CALL(c_row_anyNA, 2),    CALL(c_col_anyNA, 1),

	CALL(c_rowCollapse_init, 2),
	CALL(c_rowCollapse, 4),  CALL(c_colCollapse, 4),

	{ NULL, NULL, 0 }
};


void attribute_visible R_init_SCArray(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}
