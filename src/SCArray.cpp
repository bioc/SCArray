// ===========================================================
//
// SCArray.cpp: the C++ codes for the SCArray package
//
// Copyright (C) 2020    Xiuwen Zheng
//
// This file is part of SeqArray.
//
// SeqArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SeqArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SeqArray.
// If not, see <http://www.gnu.org/licenses/>.

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_GDS.h>

#include <set>
#include <algorithm>



// ===========================================================
// Library Functions
// ===========================================================

extern "C"
{

// ===========================================================
// the initial function when the package is loaded
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SC_ExternalName0()
{
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP SC_ExternalName1(SEXP x)
{
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP SC_ExternalName2(SEXP x, SEXP y)
{
	return R_NilValue;
}


/// initialize the package
COREARRAY_DLL_EXPORT void R_init_SCArray(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	static R_CallMethodDef callMethods[] =
	{
		CALL(SC_ExternalName0, 0),         CALL(SC_ExternalName1, 1),
		CALL(SC_ExternalName2, 2),
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Init_GDS_Routines();
}

} // extern "C"
