// ===============================================================
//
// SCArray R package
// Copyright (C) 2021   Xiuwen Zheng (@AbbVie-ComputationalGenomics)
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

SEXP c_rowVars_update(SEXP mat, SEXP pm)
{
	// get nrow and ncol
	SEXP dm = GET_DIM(mat);
	if (Rf_isNull(dm) || Rf_length(dm) != 2)
		Rf_error("It should be a numeric matrix.");
	int *pdm = INTEGER(dm);
	int nrow=pdm[0], ncol=pdm[1];
	// do
	double v, *po = REAL(pm);
	double *pS=po, *pS2=po+nrow, *pN=po+2*nrow;
	const double *p = REAL(mat);
	for (int i=0; i < ncol; i++, p+=nrow)
	{
		for (int j=0; j < nrow; j++)
		{
			if (R_FINITE(v=p[j]))
			{
				pS[j] += v; pS2[j] += v*v; 
				pN[j] ++;
			}
		}
	}
	// output
	return pm;
}

SEXP c_rowVars_final(SEXP pm)
{
	const size_t n = Rf_xlength(pm) / 3;
	const double *po = REAL(pm);
	const double *pS=po, *pS2=po+n, *pN=po+2*n;
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	for (size_t i=0; i < n; i++)
		p[i] = (pS2[i] - pS[i]*pS[i]/pN[i]) / (pN[i]-1);
	return ans;
}

}
