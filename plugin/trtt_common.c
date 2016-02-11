/*
 * TR-TT : iTeraTive Reconstruction for Tomography
 *
 * trtt_common.c --
 *
 * Basic embedded functions and tools for TR-TT.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011, the MiTiV Team.
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
 * computer knowledge. Users are ther efore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

/* PREPROCESSOR INCLUSIONS ============================================== */

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>

#include "trtt_common.h"

/* PREPROCESSOR MACROS ================================================== */

/* #define ... */

/* PRIVATE FUNCTIONS ==================================================== */

/* static ... */

/* PUBLIC FUNCTIONS ===================================================== */

int trtt_mvmult(double *out, double *in, long number, double *a, long *i, long *j)
{
  long k;
  for (k=0; k<number; ++k) {
    out[i[k]-1] += a[k]*in[j[k]-1];
  }
  return TRTT_SUCCESS;
}

/*
 * This function returns "index generator" list -- an array of longs running
 * from 1 to N, inclusive.
 */
double *_trtt_indgen(const long n)
{
  long i;
  double *x = TRTT_NEW_ARRAY(double, n);

  if (x != NULL) {
    for (i=0; i<n; ++i) {
      x[i] = (double)i;
    }
  }

  return x;
}

/*
 * This function returns the integration of array T of length N and sampling
 * STEP from bound B1 to B2, using the trapezoidal rule.
 */
double _trtt_integ(const double *t, const long n, const double step, double b1, double b2)
{
  int i;
  double result = 0;
  /* Linear interpolation factor */
  double c;
  double sc = 0.5*step;

  b1 = TRTT_MAX(0,b1);
  b2 = TRTT_MIN(n-1,b2);
  
  long b1_ceil = ceil(b1);
  long b2_floor = floor(b2);

  /* Integrate the non-truncated trapezoids */
  if (b1_ceil != b2_floor) {
    for (i=b1_ceil; i<b2_floor; ++i)
    {
      result += sc*(t[i]+t[i+1]);
    }   
  }

  /* Integrate the "left side" bound (needs a linear interpolation) */
  if (b1 != b1_ceil) {
    c = TRTT_ABS(b1_ceil - b1);
    result += sc*(c*t[b1_ceil-1] + (2-c)*t[b1_ceil]); 
  }

  /* Integrate the "right side" bound (needs a linear interpolation) */
  if (b2 != b2_floor) {
    c = TRTT_ABS(b2 - b2_floor);
    result += sc*(c*t[b2_floor+1] + (2-c)*t[b2_floor]);
  }

  return result;
}

/*
 * This function returns the partial sums of the array Y of size M and
 * constant step STEP and store them in YCUM.
 */
void _trtt_cumsum(double *ycum, const double *y, const long m, const double step)
{
  int i;
  ycum[0] = 0.0;
  for (i=1; i<m; ++i) {
    ycum[i] = ycum[i-1]+0.5*step*(y[i-1]+y[i]);
    /* printf("%f\n",ycum[i]); */
  }
}

/*
 * This function performs the same integration algorithm as Munro's one for
 * the INTEG function in Yorick. Y is the array of M points to integrate from
 * its initial point to NP regularly spaced bounds, starting from XP_I_START
 * (expressed in index units of Y), with an increment INC between points. YCUM
 * is a M-points array giving the partial sums of Y, considering a sampling
 * STEP between points.
 */
void _trtt_integ_munro(double *yp, const double *y, const double *ycum, const double xp_i_start, const double inc, const long np, const long m, const double step)
{
  double c0, c1, dx, xp_i;
  long ip, i;
  
  /* Compute interpolated integrals */
  for (i=0; i<np; ++i) {
    xp_i = xp_i_start + i*inc;
    ip = (long)floor(xp_i);
    if (ip>=0 && ip<m-1) {
      dx = xp_i-ip;
      c1 = 0.5*dx*dx;
      c0 = dx-c1;
    }
    
    if (ip<0) {
      yp[i] = 0.0;
    } else if (ip<m-1) {
      yp[i] = ycum[ip];
      yp[i] += step*(c0*y[ip]+c1*y[ip+1]);
    } else {
      yp[i] = ycum[m-1];
    } 
  }

  /* fprintf(stderr, "yp[np-1] = %f\n",yp[np-1]); */
}

/*
 * This function returns the minimum value of an array of doubles.
 */
double _trtt_get_min(const double x[], const long n)
{
  int i;
  double xmin = x[0];

  for (i=1; i<n; ++i) {
    if (x[i]<xmin) xmin = x[i];
  }

  return xmin;
}

/*
 * This function returns the maximum value of an array of doubles.
 */
double _trtt_get_max(const double x[], const long n)
{
  int i;
  double xmax = x[0];

  for (i=1; i<n; ++i) {
    if (x[i]>xmax) xmax = x[i];
  }

  return xmax;
}

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
