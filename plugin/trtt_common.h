/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * trtt_common.h --
 *
 * General macros definitions.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011, MiTiV Team.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

/* PREPROCESSOR =============================================================== */

#ifndef TRTT_COMMON_H
#define TRTT_COMMON_H 1

/* PREPROCESSOR INCLUSIONS ==================================================== */

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>

/* PREPROCESSOR MACROS ======================================================== */

/* _TRTT_BEGIN_DECLS should be used at the beginning of your declarations,
   so that C++ compilers don't mangle their names.  Use _TOMO_END_DECLS at
   the end of C declarations. */

/*
 * Title: TRTT common macros and functions
 *
 *   Version 1.0
 *
 *   Copyright (C) 2011, the MiTiV project.
 *
 *   This package defines macros and functions commonly used in TR-TT
 *   packages.
 *   The definitions are available by: 
 *   > #include <trtt_common.h>
 */

/*
 * Defines: Miscellaneous macros.
 *
 *   TRTT_BEGIN_C_DECLS - Start declarations of C elements (functions, variables,
 *                        types, structures, etc.).
 *   TRTT_END_C_DECLS   - Terminate C declarations.
 *
 *   <TRTT_BEGIN_C_DECLS> should be used at the beginning of your declarations,
 *   so that C++ compilers don't mangle their names.  Use <TRTT_END_C_DECLS> at
 *   the end of C declarations.
 */
#ifdef __cplusplus
# define TRTT_BEGIN_C_DECLS extern "C" {
# define TRTT_END_C_DECLS }
#else /* not __cplusplus */
# define TRTT_BEGIN_C_DECLS /* empty */
# define TRTT_END_C_DECLS /* empty */
#endif /* __cplusplus */

/*
 * Macro: TRTT_TRUE
 *   Boolean value symbolizing a TRUE answer to a test.
 *
 * Description: 
 *   <TRTT_TRUE> and <TRTT_FALSE> are the standard values used to answer YES
 *   or NO to a condition.
 *
 * See Also:
 *   <TRTT_FALSE>.
 */
#define TRTT_TRUE 1

/*
 * Macro: TRTT_FALSE
 *   Boolean value symbolizing a FALSE answer to a test.
 *
 * Description: 
 *   <TRTT_TRUE> and <TRTT_FALSE> are the standard values used to answer YES
 *   or NO to a condition.
 *
 * See Also:
 *   <TRTT_TRUE>.
 */
#define TRTT_FALSE 0

/*
 * Macro: TRTT_SUCCESS
 *   Returned value upon success.
 *
 * Description:
 *   <TRTT_SUCCESS> and <TRTT_FAILURE> are the standard values returned by
 *   most dynamic buffer routines.
 *
 * See Also:
 *   <TRTT_FAILURE>.
 */
#define TRTT_SUCCESS   (0)

/*
 * Macro: TRTT_FAILURE
 *   Returned value upon failure.
 *
 * See Also:
 *   <TRTT_SUCCESS>.
 */
#define TRTT_FAILURE   (-1)

/*
 * Macro: TRTT_WRONG_VALUE 
 *   Returned value upon failure for functions returning a scalar value.
 *
 */
#define TRTT_WRONG_VALUE   (-1)

/*
 * Macro: TRTT_FOOTPRINT_SIZE 
 * Number of points for the storage of the normalized footprint on
 * detector over one dimension.
 */
#define TRTT_SIZE_FOOTPRINT (1000)

/*
 * Macro: TRTT_SPLINE_DEG
 * Default degree of the spatial Bsplines for the object model.
 */
#define TRTT_SPLINE_DEG (3)

/*
 * Macro: TRTT_T_SPLINE_DEG
 * Default degree of the spatial Bsplines for the object model.
 */
#define TRTT_T_SPLINE_DEG (3)

/*
 * Macro: TRTT_NEW
 *
 * TRTT_NEW(_type_) allocates an object of structure/class _type_.
 *
 * See Also:
 *   <TRTT_NEW0>, <TRTT_NEW_ARRAY>.
 */
#define TRTT_NEW(type)               ((type *)malloc(sizeof(type)))

/*
 * Macro: TRTT_NEW0
 *
 * TRTT_NEW0(_type_) allocates an object of structure/class _type_, and fills
 * it with zero values.
 *
 * See Also:
 *   <TRTT_NEW>, <TRTT_NEW_ARRAY>.
 */
#define TRTT_NEW0(type)              ((type *)calloc(1, sizeof(type)))

/*
 * Macro: TRTT_NEW_ARRAY
 *
 * TRTT_NEW_ARRAY(_type_) allocates an array of _number_ items of
 * structure/class _type_.
 *
 * See Also:
 *   <TRTT_NEW>, <TRTT_NEW0>.
 */
#define TRTT_NEW_ARRAY(type, number) ((type *)malloc((number)*sizeof(type)))

/*
 * Macro: TRTT_NEW_ARRAY0
 *
 * TRTT_NEW_ARRAY0(_type_) allocates an array of _number_ items of
 * structure/class _type_, and fills it with zero values.
 *
 * See Also:
 *   <TRTT_NEW>, <TRTT_NEW0>.
 */
#define TRTT_NEW_ARRAY0(type, number) ((type *)calloc((number),sizeof(type)))

/*
 * Macro: TRTT_ABS
 *
 * TRTT_ABS(_a_) gives the absolute value of _a_.
 */
#define TRTT_ABS(a)   ((a) >= 0 ? (a) : -(a))

/*
 * Macro: TRTT_MIN
 *
 * TRTT_MIN(_a_, _b_) yields the minimum value between _a_ and _b_.
 */
#define TRTT_MIN(a,b)                ((a) <= (b) ? (a) : (b))

/*
 * Macro: TRTT_MAX
 *
 * TRTT_MAX(_a_, _b_) yields the maximum value between _a_ and _b_.
 */
#define TRTT_MAX(a,b)                ((a) >= (b) ? (a) : (b))

/*
 * Macro: TRTT_PIXEL_VALUE_ACCESS
 *
 * TRTT_PIXEL_VALUE_ACCESS(_a_, _nx_, _i_, _j_) gives the value of the pixel
 * situated at coordinates (_i_, _j_) in the 2D image _a_ assimilated to a
 * matrix _ny_ rows x _nx_ columns (row major storage).
 */
#define TRTT_PIXEL_VALUE_ACCESS(a,nx,i,j)        ((a)[(j)*(nx) + (i)])

/*
 * Macro: TRTT_VOXEL_VALUE_ACCESS
 *
 * TRTT_VOXEL_VALUE_ACCESS(_a_, _nx_, _ny_, _i_, _j_, _k_) gives the value of the voxel
 * situated at coordinates (_i_, _j_, _k_) in the 3D image _a_ assimilated to a
 * matrix _nz_ x _ny_ x _nx_ (row major storage).
 */
#define TRTT_VOXEL_VALUE_ACCESS(a,nx,ny,i,j,k)   ((a)[((k)*(ny) +(j))*(nx) + (i)])

/*
 * Macro: TRTT_T_VOXEL_VALUE_ACCESS
 *
 * TRTT_VOXEL_VALUE_ACCESS(_a_, _nx_, _ny_, _nz_, _i_, _j_, _k_, _l_) gives the value
 * of the spatio-temporal voxel situated at coordinates (_i_, _j_, _k_, _l_) in
 * the 4D image _a_ assimilated to a matrix _nt_ x _nz_ x _ny_ x _nx_ (row
 * major storage).
 */
#define TRTT_T_VOXEL_VALUE_ACCESS(a,nx,ny,nz,i,j,k,l)   ((a)[(nx)*((ny)*((l)*(nz) + (k)) + (j)) + (i)])

TRTT_BEGIN_C_DECLS

/* PROTOYPES ================================================================== */

/*
 * This function applies a sparse matrix multiplication.
 */
extern int trtt_mvmult(double *out, double *in, long number, double *a, long *i, long *j);

/*
 * This function returns "index generator" list -- an array of longs running
 * from 1 to N, inclusive.
 */
extern double *_trtt_indgen(const long n);

/*
 * This function returns the integration of array T of length N and sampling
 * STEP from bound B1 to B2, using the trapezoidal rule.
 */
extern double _trtt_integ(const double *t, const long n, const double step, double b1, double b2);

/*
 * This function returns the partial sums of the array Y of size M and
 * constant step STEP and store them in YCUM.
 */
extern void _trtt_cumsum(double *ycum, const double *y, const long m, const double step);

/*
 * This function performs the same integration algorithm as Munro's one for
 * the INTEG function in Yorick. Y is the array of M points to integrate from
 * its initial point to NP regularly spaced bounds, starting from XP_I_START
 * (expressed in index units of Y), with an increment INC between points. YCUM
 * is a M-points array giving the partial sums of Y, considering a sampling
 * STEP between points.
 */
extern void _trtt_integ_munro(double *yp, const double *y, const double *ycum, const double xp_i_start, const double inc, const long np, const long m, const double step);

/*
 * This function returns the minimum value of an array of doubles.
 */
extern double _trtt_get_min(const double x[], const long n);

/*
 * This function returns the maximum value of an array of doubles.
 */
extern double _trtt_get_max(const double x[], const long n);

TRTT_END_C_DECLS

#endif /* TRTT_COMMON_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 78
 * coding: utf-8
 * End:
 */
