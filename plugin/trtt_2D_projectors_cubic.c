/*
 * TR-TT : iTeraTive Reconstruction for Tomography
 * 
 *
 * trtt_2D_projectors.c --
 *
 * TRTT 2D projectors.
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

/* PREPROCESSOR INCLUSIONS ============================================== */

#define ONE_TWENTY_FOURTH 0.0416666666666667
#define ONE_THIRD 0.3333333333333333
#define TWO_THIRD 0.6666666666666666
#define FOUR_THIRD 1.3333333333333333
#define ONE_EIGHTH 0.125

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "trtt_common.h"
#include "trtt_xform3d.h"
#include "trtt_2D_projectors.h"


/* PREPROCESSOR MACROS ================================================== */

#define STATEMENT(code) do { code; } while (0)

/* PRIVATE FUNCTIONS ==================================================== */

////////////////////////////////////////
// PRIVATE FUNCTIONS FOR 2D OPERATORS //
////////////////////////////////////////

static double _trtt_cubic_bspline_integrale(double x)
{
  double x_sqr = x*x;
  double x_cub = x_sqr*x;
  double y;

  if (x<=-2.0) {
    y = 0.0;
  } else if (x>-2.0 && x<=-1.0) {
    y = ONE_TWENTY_FOURTH*(x*x_cub) + ONE_THIRD*x_cub + x_sqr + FOUR_THIRD*x + TWO_THIRD;
  } else if (x>-1.0 && x<=0.0) {
    y = (-1.0)*ONE_EIGHTH*(x*x_cub) - ONE_THIRD*x_cub + TWO_THIRD*x + 0.5;
  } else if (x>0.0 && x<=1.0) {
    y = ONE_EIGHTH*(x*x_cub) - ONE_THIRD*x_cub + TWO_THIRD*x + 0.5;
  } else if (x>1.0 && x<2.0) {
    y = (-1.0)*ONE_TWENTY_FOURTH*(x*x_cub) + ONE_THIRD*x_cub - x_sqr + FOUR_THIRD*x + ONE_THIRD;
  } else {
    y = 1.0;
  }

  return y;
}

/*----------------*/
/*-Burn footprint-*/
/*----------------*/
static int _trtt_2D_burn_footprint(double *pix, 
				   const double vox,
				   const double vc,
				   const double v_left,
				   const double v_right,
				   const double rho0,
				   const double h_norm,
				   const double scl_norm)
{
  long i;
  double p_first, v_next, p_next;
  long npix;
  /* double coord_footprint_0; */
  double coord_footprint_start;
  
  long v_left_id = (long)v_left;
  // FIXME: we are in bound detector pixel coordinates
  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv
  
  if (v_left<v_right) {
    /* Get the position of left-hand side of the footprint in pixel detector
     * space */
    coord_footprint_start = (v_left-vc)*scl_norm;
    // FIXME: generic footprint coordinates for cubic b-spline
    // FIXME: -2 ; -1 ; 0 ; 1 ; 2
    
    npix = (long)(v_right-v_left); // FIXME: number of impacted detector pixels
    /*FIXME: CUBIC B-SPLINE OPTIMIZED CALCULATION */
    p_first = _trtt_cubic_bspline_integrale(coord_footprint_start);
    
    for (i=0; i<npix; ++i) {    
      v_next = coord_footprint_start + (i+1)*scl_norm;
      p_next = _trtt_cubic_bspline_integrale(v_next);
      
      /* Accumulate in current detector pixel */
      *(pix+v_left_id+i) += (p_next-p_first)*vox*rho0;
      //FIXME: EXECUTION TIME BOTTLENECK !!!
      
      p_first = p_next;
    }
  }

  return TRTT_SUCCESS;
}
/*----------------------------*/
/*-Burn footprint (transpose)-*/
/*----------------------------*/
static int _trtt_2D_burn_footprint_transp(double *vox,
					  const double *pix,  
					  const double vc,  
					  const double v_left,
					  const double v_right,
					  const long nx,
					  const long id_x,
					  const long id_y,
					  const double rho0,
					  const double h_norm,
					  const double scl_norm)
{
  long i;
  double p_first, v_next, p_next;
  long npix;
  /* double coord_footprint_0; */
  double coord_footprint_start;
  /* Get the pointer to the current voxel */
  double *current_voxel = vox + nx*id_y + id_x;
  long v_left_id = (long)v_left;
  // FIXME: we are in bound detector pixel coordinates
  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv
  
  if (v_left < v_right) {
    /* Get the position of left-hand side of the footprint in pixel detector
     * space */
    coord_footprint_start = (v_left-vc)*scl_norm;
    // FIXME: generic footprint coordinates for cubic b-spline
    // FIXME: -2 ; -1 ; 0 ; 1 ; 2
    
    npix = (long)(v_right-v_left); // FIXME: number of impacted detector pixels
    
    /*FIXME: CUBIC B-SPLINE OPTIMIZED CALCULATION */
    p_first = _trtt_cubic_bspline_integrale(coord_footprint_start);
    
    for (i=0; i<npix; ++i) {
      v_next = coord_footprint_start + (i+1)*scl_norm;
      p_next = _trtt_cubic_bspline_integrale(v_next);
      
      /* Accumulate in current detector pixel */
      *current_voxel += rho0 * (*(pix+v_left_id+i)) * (p_next - p_first);
      //FIXME: EXECUTION TIME BOTTLENECK !!!
    
      p_first = p_next;
    }
  }

  return TRTT_SUCCESS;
}

/* PUBLIC FUNCTIONS =================================================== */

/*---------------------------------------------------------------*/
/* See NaturalDocs Documentation for the help of these functions */
/*---------------------------------------------------------------*/

//////////////////
// 2D OPERATORS //
//////////////////

/*------------------*/
/*-2D PARALLEL BEAM-*/
/*------------------*/
int trtt_PB2D_cubic(double *out, double *in, long nx, long ny, double s_scl, int s_deg, long nv, double v_scl, double *xos, double *xsd, int job)
{
  /*  Calculated parameters */
  double in_ij;
  double vc, vc_temp, v_left, v_right;
  long i,j;
  /* Extract some XSD coefficients */
  double a1 = xsd[0];
  double a1inv = -1.0*(1.0/a1);
  double a2 = xsd[1];
  double a4 = xsd[3];
  double a5 = xsd[4];
  double a5a1 = a5*a1inv; 
  double a5a1scl = a5a1*v_scl; //FIXME: metric scale
  double a6 = xsd[5];
  double a8 = xsd[7];
  double b1 = a6 + a5a1*a2;
  double b3 = a8 + a5a1*a4;
  double b1o5 = b1*xos[4];
  double b1o6 = b1*xos[5];
  double b1o8b3 = b1*xos[7] + b3;

  /* Calculate the distorsion factor */
  double delta_v = sqrt(1+a5a1scl*a5a1scl);
  double rho0 = s_scl*s_scl*delta_v;
  /* Normalized b-spline step to pixel detector sampling  */
  double h_norm = (delta_v*s_scl)/v_scl;
  double scl_norm = 1/h_norm; //FIXME: scl_norm = 1/(h_norm*step_footprint);
  /* Retrieve the half support of the footprint on detector */
  double half_support = 0.5*(s_deg+1)*h_norm;

  for (j=0; j<ny; ++j) {
      /* Pre-projection on Detector */
      vc_temp = b1o6*j + b1o8b3;
      
      for (i=0; i<nx; ++i) {	  
	  /* Projection on Detector */
	  vc = vc_temp + b1o5*i;
	  
	  /* Get the bounds positions of the footprint */
          v_left = vc-half_support;
          v_right = vc+half_support;

          if (v_left < 0.0) v_left = 0.0;
          else if (v_left >= 0.0 && v_left <= (nv-1)) v_left = floor(v_left);
          else v_left = nv-1;

          if (v_right < 0.0) v_right = 0.0;
          else if (v_right >= 0.0 && v_right <= (nv-2)) v_right = floor(v_right)+1;
          else v_right = nv-1;
	  /* v_left = floor(vc-half_support); */
	  /* v_right = floor(vc+half_support); */
	  /* v_left = TRTT_MIN(TRTT_MAX(v_left,0),nv-1); */
	  /* v_right = TRTT_MAX(TRTT_MIN(v_right,nv-1),0); */
	  // FIXME: we are in bound detector pixel coordinates
	  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv
	  
	  /* Burn footprint on Detector */
	  if (job == 0) {
	      /* DIRECT OPERATOR */
	      /* Get voxel value */
	      in_ij = TRTT_PIXEL_VALUE_ACCESS(in,nx,i,j);
              if (_trtt_2D_burn_footprint(out, in_ij, vc, v_left, v_right, \
                                          rho0, h_norm, scl_norm) == TRTT_FAILURE) {
                return TRTT_FAILURE;
              }   
	  } else {
	      /* TRANSPOSE OPERATOR */
	      if (_trtt_2D_burn_footprint_transp(out, in, vc, v_left, v_right, 
						 nx, i, j,		\
						 rho0, h_norm, scl_norm) == TRTT_FAILURE) {
		  return TRTT_FAILURE;
	      }
	  }
      }
  }
  
  return TRTT_SUCCESS;
}

/*-------------*/
/*-2D FAN BEAM-*/
/*-------------*/
int trtt_FB2D_cubic(double *out, double *in, long nx, long ny, double s_scl, int s_deg, long nv, double v_scl, double *xos, double *xsd, int job)
{
    /*  Calculated parameters */
    double in_ij;
    double wk, wk_temp, vc, v_left, v_right, vca8;
    double num, den, num_temp, den_temp;
    long i,j;
    /* Extract some XSD coefficients */
    double a1 = xsd[0];
    double a1inv = -1.0*(1.0/a1);
    double a2 = xsd[1];
    double a4 = xsd[3];
    double a4m = -1.0*a4;
    double a4sqr = a4*a4;
    double a5 = xsd[4];
    double a6 = xsd[5];
    double a8 = xsd[7];
    double b1 = a6 + a1inv*a5*a2;
    double b3 = a8 + a1inv*a5*a4;

    double b1o5 = b1*xos[4];
    double b1o6 = b1*xos[5];
    double b1o8 = b1*xos[7];
    double a1o1a2o5 = a1*xos[0] + a2*xos[4];
    double a1o2a2o6 = a1*xos[1] + a2*xos[5];
    double a1o4a2o8 = a1*xos[3] + a2*xos[7];

    /* Calculate the distorsion factor */
    double magn;
    double delta_v;
    double rho0;
    double rho0_init = s_scl*s_scl;
    /* Normalized b-spline step to pixel detector sampling  */
    double h_norm;
    double h_norm_init = s_scl/v_scl;
    double scl_norm;
    /* Retrieve the half support of the footprint on detector */
    double half_support;
    double half_support_init = 0.5*(s_deg+1);

    for (j=0; j<ny; ++j) {
        /* Pre-projection on Detector */
        num_temp = b1o6*j + b1o8;
        den_temp = a1o2a2o6*j + a1o4a2o8;
        wk_temp = xos[1]*j + xos[3];

	for (i=0; i<nx; ++i) {
	    /* Projection on Detector */
	    wk = wk_temp + xos[0]*i;
            num = num_temp + b1o5*i;
            den = den_temp + a1o1a2o5*i;
            vc = (a4m/den)*num + b3;
            
	    /* Calculate the distorsion factors */    
	    if (wk != 0) { 
              magn = a4m/wk;         
              vca8 = (vc-a8)*v_scl; //FIXME: metric scale
              delta_v = sqrt(1.0+((vca8*vca8)/(a4sqr)));
              rho0 = rho0_init * magn * delta_v;
              h_norm = h_norm_init*magn*delta_v;
              scl_norm = 1.0/h_norm;//FIXME: scl_norm = 1/(h_norm*step_footprint);
              half_support = half_support_init * h_norm;
	    } else {
              fprintf(stderr, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n");
              return TRTT_FAILURE;
	    }

	    /* Get the bounds positions of the footprint */
            v_left = vc-half_support;
            v_right = vc+half_support;

            if (v_left < 0.0) v_left = 0.0;
            else if (v_left >= 0.0 && v_left <= (nv-1)) v_left = floor(v_left);
            else v_left = nv-1;
            
            if (v_right < 0.0) v_right = 0.0;
            else if (v_right >= 0.0 && v_right <= (nv-2)) v_right = floor(v_right)+1;
            else v_right = nv-1;
	    /* v_left = floor(vc-half_support); */
	    /* v_right = floor(vc+half_support); */
	    /* v_left = TRTT_MIN(TRTT_MAX(v_left,0),nv-1); */
	    /* v_right = TRTT_MAX(TRTT_MIN(v_right,nv-1),0); */
	    // FIXME: we are in bound detector pixel coordinates
	    // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv
	    
	    /* Burn footprint on Detector */
	    if (job == 0) {
              /* DIRECT OPERATOR */
              /* Get voxel value */
              in_ij = TRTT_PIXEL_VALUE_ACCESS(in,nx,i,j);
              if (_trtt_2D_burn_footprint(out, in_ij, vc, v_left, v_right, \
                                          rho0, h_norm, scl_norm) == TRTT_FAILURE) {
                return TRTT_FAILURE;
              }  
	    } else {
		/* TRANSPOSE OPERATOR */
		if (_trtt_2D_burn_footprint_transp(out, in, vc, v_left, v_right, 
						   nx, i, j,	\
						   rho0, h_norm, scl_norm) == TRTT_FAILURE) {
		    return TRTT_FAILURE;
		}
	    }
	}
    }
    
    return TRTT_SUCCESS;
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
