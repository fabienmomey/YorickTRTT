/*
 * TR-TT : iTeraTive Reconstruction for Tomography
 * 
 *
 * trtt_projectors.c --
 *
 * TRTT Radon projectors.
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

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "trtt_common.h"
#include "trtt_xform3d.h"
#include "trtt_projectors.h"


/* PREPROCESSOR MACROS ================================================== */

#define STATEMENT(code) do { code; } while (0)

/* PRIVATE FUNCTIONS ==================================================== */

////////////////////////////////////////
// PRIVATE FUNCTIONS FOR 2D OPERATORS //
////////////////////////////////////////

/*----------------------------*/
/*-Object to Source Transform-*/
/*----------------------------*/
static int _trtt_2D_transform_O_to_S(double *coord, const double transform_coefs[12])
{
  double outcoord[2];
  double x = coord[0];
  double y = coord[1];

  /* Calculate output coordinates */
  outcoord[0] = transform_coefs[0]*x + transform_coefs[1]*y + transform_coefs[3];
  outcoord[1] = transform_coefs[4]*x + transform_coefs[5]*y + transform_coefs[7];

  memcpy(coord, outcoord, 2*sizeof(double));

  return TRTT_SUCCESS;
}
/*---------------------------------*/
/*-Parallel projection on Detector-*/
/*---------------------------------*/
static int _trtt_2D_projection_S_to_D_parallel(double *coord, const double transform_coefs[12])
{
  double outcoord[2];
  /* double w0 = coord[0]; */
  double v0 = coord[1];

  double a1 = transform_coefs[0];  
  if (a1 == 0) {
    fprintf(stderr, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n");
    return TRTT_FAILURE;
  }
  
  double a1inv = -1.0*(1.0/a1);
  double a2 = transform_coefs[1];
  double a4 = transform_coefs[3];
  double a5 = transform_coefs[4];
  double a6 = transform_coefs[5];
  double a8 = transform_coefs[7];

  double b1 = a6 + a1inv*a5*a2;
  double b3 = a8 + a1inv*a5*a4;

  outcoord[0] = a1inv*(a2*v0 + a4);
  outcoord[1] = b1*v0 + b3;

  memcpy(coord, outcoord, 2*sizeof(double));

  return TRTT_SUCCESS;
}
/*----------------------------*/
/*-Fan projection on Detector-*/
/*----------------------------*/
static int _trtt_2D_projection_S_to_D_fan(double *coord, const double transform_coefs[12], double tan_gamma0)
{
  double outcoord[2];
  /* double v0 = coord[1]; */

  double a1 = transform_coefs[0];
  if (a1 == 0) {
    fprintf(stderr, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n");
    return TRTT_FAILURE;
  }
  
  double a1inv = -1.0*(1.0/a1);
  double a2 = transform_coefs[1];
  double a4 = transform_coefs[3];
  double a5 = transform_coefs[4];
  double a6 = transform_coefs[5];
  double a8 = transform_coefs[7];

  double b1 = a6 + a1inv*a5*a2;
  double b3 = a8 + a1inv*a5*a4;

  double w = -(1.0/(a1+a2*tan_gamma0))*a4;
  outcoord[0] = w;
  outcoord[1] = b1*tan_gamma0*w + b3;

  memcpy(coord, outcoord, 2*sizeof(double));

  return TRTT_SUCCESS;
}
/*----------------*/
/*-Burn footprint-*/
/*----------------*/
static int _trtt_2D_burn_footprint(double *pix, 
                                   const double vox,
                                   const double *coord,
                                   const double v_scl,
                                   const long nv,
                                   const double s_scl, 
                                   const int s_deg,
                                   const double rho_0,
                                   const double *footprint, 
                                   const double *cum_footprint, 
                                   const long size_footprint,
                                   const double step_footprint,
                                   const double *coord_footprint)
{
  long i;
  double *pix_integ;
  long npix;
  double coord_footprint_0;
  double coord_footprint_start;
  /* Normalized b-spline step to pixel detector sampling  */
  double h_norm = (1.0/v_scl)*(rho_0/s_scl);
  double scl_norm = 1/(h_norm*step_footprint);
  /* Retrieve the half support of the footprint on detector */
  double half_support = 0.5 * (s_deg+1) * h_norm;
  /* Get the position of impact point of the blob center */
  double vc =  coord[1];
  /* Get the bounds positions of the footprint */
  double v_left = vc-half_support;
  double v_right = vc+half_support;
  // FIXME: we are in bound detector pixel coordinates
  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv

  /* Get the indexes of impacted pixels detector */
  double d_left = floor(v_left);
  double d_right = floor(v_right);
  d_left = TRTT_MIN(TRTT_MAX(d_left,0),nv-1);
  d_right = TRTT_MAX(TRTT_MIN(d_right,nv-1),0);
  
  if (d_left <= d_right) {
    /* Get the position of left-hand side of the footprint in pixel detector
     * space */
    coord_footprint_0 = vc + coord_footprint[0]*h_norm;
    coord_footprint_start = (d_left-coord_footprint_0)*scl_norm;
    // FIXME: footprint samples coordinates
    // FIXME: 0 ; 1 ; 2 ; ... ; size_footprint-1
    
    npix = d_right-d_left+1; // FIXME: number of impacted detector pixels
    /* fprintf(stderr, "d_left = %i    d_right = %i\n",(int)d_left,(int)d_right); */

    // FIXME: bound effects are handled by INTEG function

    pix_integ = TRTT_NEW_ARRAY0(double, npix+1);
    if (pix_integ == NULL) {
      return TRTT_FAILURE;
    }

    _trtt_integ_munro(pix_integ, footprint,                             \
                      cum_footprint, coord_footprint_start,             \
                      scl_norm, npix+1, size_footprint, step_footprint);

    for (i=0; i<npix; ++i) {
      *(pix+(long)d_left+i) += rho_0*vox*(pix_integ[i+1]-pix_integ[i]);
    }

    free(pix_integ);
  }

  return TRTT_SUCCESS;
}
/*----------------------------*/
/*-Burn footprint (transpose)-*/
/*----------------------------*/
static int _trtt_2D_burn_footprint_transp(double *vox,
                                          const double *pix,    
                                          const double *coord,  
                                          const double v_scl,
                                          const long nv,
                                          const long nx,
                                          const long ny,
                                          const long id_x,
                                          const long id_y,
                                          const double s_scl, 
                                          const int s_deg,
                                          const double rho_0,
                                          const double *footprint, 
                                          const double *cum_footprint,  
                                          const long size_footprint,
                                          const double step_footprint,
                                          const double *coord_footprint)
{
  long i;
  double *pix_integ;
  long npix;
  double coord_footprint_0;
  double coord_footprint_start;
  /* Get the pointer to the current voxel */
  double *current_voxel = vox + nx*id_y + id_x;
  /* Normalized b-spline step to pixel detector sampling  */
  double h_norm = (1.0/v_scl)*(rho_0/s_scl);
  double scl_norm = 1.0/(step_footprint*h_norm);
  /* Retrieve the half support of the footprint on detector */
  double half_support = 0.5 * (s_deg+1) * h_norm;
  /* Get the position of impact point of the blob center */
  double vc =  coord[1];
  /* Get the bounds positions of the footprint */
  double v_left = vc-half_support;
  double v_right = vc+half_support;

  /* Get the indexes of impacted pixels detector */
  double d_left = floor(v_left);
  double d_right = floor(v_right);
  d_left = TRTT_MIN(TRTT_MAX(d_left,0),nv-1);
  d_right = TRTT_MAX(TRTT_MIN(d_right,nv-1),0);
  // FIXME: we are in bound detector pixel coordinates
  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv

  if (d_left <= d_right) {
    /* Get the position of left-hand side of the footprint in pixel detector
     * space */
    coord_footprint_0 = vc + coord_footprint[0]*h_norm;
    coord_footprint_start = ((double)d_left-coord_footprint_0)*scl_norm;
    // FIXME: footprint samples coordinates
    // FIXME: 0 ; 1 ; 2 ; ... ; size_footprint-1

    npix = d_right-d_left+1; // FIXME: number of impacted detector pixels

    // FIXME: bound effects are handled by INTEG function

    pix_integ = TRTT_NEW_ARRAY0(double, npix+1);
    if (pix_integ == NULL) {
      return TRTT_FAILURE;
    }

    _trtt_integ_munro(pix_integ, footprint,                             \
                      cum_footprint, coord_footprint_start,             \
                      scl_norm, npix+1, size_footprint, step_footprint);
      
    
    for (i=0; i<npix; ++i) {
      *current_voxel += rho_0 * (*(pix+(long)d_left+i)) * (pix_integ[i+1]-pix_integ[i]);
    }

    free(pix_integ);
  }

  return TRTT_SUCCESS;
}

////////////////////////////////////////
// PRIVATE FUNCTIONS FOR 3D OPERATORS //
////////////////////////////////////////

/* static int _trtt_3D_transform_O_to_S(double *coord, const double transform_coefs[12]) */
/* { */
/*   return TRTT_FAILURE; */
/* } */

/* static int _trtt_3D_projection_S_to_D_parallel(double *coord, const double transform_coefs[12]) */
/* { */
/*   return TRTT_FAILURE; */
/* } */

/* static int _trtt_3D_projection_S_to_D_cone(double *coord, const double transform_coefs[12], double tan_gamma0, double tan_alpha, double *magn, double *delta_u, double *delta_v) */
/* { */
/*   return TRTT_FAILURE; */
/* } */

/* static int _trtt_3D_burn_footprint(double *coord, double magn, double delta_u, double delta_v) */
/* { */
/*   return TRTT_FAILURE; */
/* } */

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
int trtt_PB2D(double *out, double *in, long nx, long ny, double s_scl, double rho_0, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nv, double v_scl, double *xos, double *xsd, int job)
{
  /*  Calculated parameters */
  double in_ij;
  double coord[2];
  int i,j;

  for (j=0; j<ny; ++j) {
    for (i=0; i<nx; ++i) {
      coord[0] = i;
      coord[1] = j;
      
      /* Transform Object -> Source */
      _trtt_2D_transform_O_to_S(coord, xos);
      
      /* Projection on detector */
      if (_trtt_2D_projection_S_to_D_parallel(coord, xsd) == TRTT_FAILURE) {
        return TRTT_FAILURE;
      }
      
      if (job == 0) {
        /* DIRECT OPERATOR */
        /* Get voxel value */
        in_ij = TRTT_PIXEL_VALUE_ACCESS(in,nx,i,j);
        if (_trtt_2D_burn_footprint(out, in_ij, coord, \
                                    v_scl, nv, s_scl,                   \
                                    s_deg, rho_0,        \
                                    footprint, cum_footprint, size_footprint, \
                                    step_footprint, coord_footprint) == TRTT_FAILURE) {
          return TRTT_FAILURE;
        }      
      } else {
        /* TRANSPOSE OPERATOR */
        if (_trtt_2D_burn_footprint_transp(out, in, coord, v_scl, nv, \
                                           nx, ny, i, j, s_scl, \
                                           s_deg, rho_0,        \
                                           footprint, cum_footprint, size_footprint, \
                                           step_footprint, coord_footprint) == TRTT_FAILURE) {
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
int trtt_FB2D(double *out, double *in, long nx, long ny, double s_scl, double rho_0, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nv, double v_scl, double *xos, double *xsd, int job)
{
  /* Calculated parameters */
  double in_ij;
  double coord[2];
  double tan_gamma0, w0, rho_0_old, focale;
  int i,j;

  rho_0_old = rho_0;
  focale = xsd[3];

  for (j=0; j<ny; ++j) {
    for (i=0; i<nx; ++i) {
      coord[0] = i;
      coord[1] = j;
      
      /* Transform Object -> Source */
      _trtt_2D_transform_O_to_S(coord, xos);
      
      /* Projection on detector */
      w0 = coord[0];
      tan_gamma0 = coord[1]/w0;

      if (_trtt_2D_projection_S_to_D_fan(coord, xsd, tan_gamma0) == TRTT_FAILURE) {   
        return TRTT_FAILURE;
      }

      if (prev_x != 0) {
        rho_0 = rho_0_old * (focale/w0) * sqrt(1.0+tan_gamma0*tan_gamma0);
      } else {
        fprintf(stderr, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n");
        return TRTT_FAILURE;
      }

      /* Burn footprint on Detector */
      if (job == 0) {
        /* DIRECT OPERATOR */
        /* Get voxel value */
        in_ij = TRTT_PIXEL_VALUE_ACCESS(in,nx,i,j);
        if (_trtt_2D_burn_footprint(out, in_ij, coord, \
                                    v_scl, nv, s_scl,               \
                                    s_deg, rho_0,               \
                                    footprint, cum_footprint, size_footprint, \
                                    step_footprint, coord_footprint) == TRTT_FAILURE) {
          return TRTT_FAILURE;
        }
      } else {
        /* TRANSPOSE OPERATOR */
        if (_trtt_2D_burn_footprint_transp(out, in, coord, v_scl, nv, \
                                           nx, ny, i, j, s_scl, \
                                           s_deg, rho_0,        \
                                           footprint, cum_footprint, size_footprint, \
                                           step_footprint, coord_footprint) == TRTT_FAILURE) {
          return TRTT_FAILURE;
        }
      }
    }
  }
  
  return TRTT_SUCCESS;
}

//////////////////
// 3D OPERATORS //
//////////////////

/* int trtt_PB3D() */
/* { */
/*   return TRTT_FAILURE; */
/* } */

/* int trtt_CB3D() */
/* { */
/*   return TRTT_FAILURE; */
/* } */

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
