/*
 * TR-TT : iTeraTive Reconstruction for Tomography
 *
 * trtt_3D_projectors.c --
 *
 * TRTT 3D projectors.
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

/* #define ONE_TWENTY_FOURTH 0.0416666666666667 */
/* #define ONE_THIRD 0.3333333333333333 */
/* #define TWO_THIRD 0.6666666666666666 */
/* #define FOUR_THIRD 1.3333333333333333 */
/* #define ONE_EIGHTH 0.125 */

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "trtt_common.h"
#include "trtt_xform3d.h"
#include "trtt_3D_projectors.h"

/*
 * Title: TRTT Radon projectors for tomography (details of implementation)
 *
 *   Version 1.0
 *
 *   Copyright (C) 2011, the MiTiV project.
 *
 *   This package contains the projectors performing the algebraic Radon
 *   transform at several dimensions and geometries.
 *
 *   The definitions are available by:
 *   >  #include <trtt_3D_projectors.c>
 */

/*---------------------------------------------------------------------------*/

/* PREPROCESSOR MACROS ================================================== */

#define STATEMENT(code) do { code; } while (0)

/* PRIVATE FUNCTIONS ==================================================== */

////////////////////////////////////////
// PRIVATE FUNCTIONS FOR 3D OPERATORS //
////////////////////////////////////////

/* Function: _trtt_cubic_bspline_integrale_3D
 *
 * Description:
 *
 *      This function calculates the integrale of monodimensional B-spline of
 *      degree 3 from minus infinity to the x position
 *
 *      >           _______
 *      >          /   | ' \
 *      >         /|   | ' |\
 *      >        / |   | ' | \
 *      >       /  |   | ' |  \
 *      >  ____/   |   | ' |   \____
 *      >                '
 *      >-Inf ---------> x
 *      >
 *      >     -2  -1   0   1   2
 *      >
 *      >      <h ><h ><h ><h >
 * 
 *      The considered B-spline is normalized : h=1. As a result x has to be
 *      passed normalized in this space, i.e. centered at the center of the
 *      B-spline and divided by the scale h (in this application the sampling
 *      step of voxels).
 *
 * Parameters: 
 *
 *      x - upper integration bound.
 *
 */
/* static double _trtt_cubic_bspline_integrale_3D(double x) */
/* { */
/*   double x_sqr = x*x; */
/*   double x_cub = x_sqr*x; */
/*   double y; */

/*   if (x<=-2.0) { */
/*     y = 0.0; */
/*   } else if (x>-2.0 && x<=-1.0) { */
/*     y = ONE_TWENTY_FOURTH*(x*x_cub) + ONE_THIRD*x_cub + x_sqr + FOUR_THIRD*x + TWO_THIRD; */
/*   } else if (x>-1.0 && x<=0.0) { */
/*     y = (-1.0)*ONE_EIGHTH*(x*x_cub) - ONE_THIRD*x_cub + TWO_THIRD*x + 0.5; */
/*   } else if (x>0.0 && x<=1.0) { */
/*     y = ONE_EIGHTH*(x*x_cub) - ONE_THIRD*x_cub + TWO_THIRD*x + 0.5; */
/*   } else if (x>1.0 && x<2.0) { */
/*     y = (-1.0)*ONE_TWENTY_FOURTH*(x*x_cub) + ONE_THIRD*x_cub - x_sqr + FOUR_THIRD*x + ONE_THIRD; */
/*   } else { */
/*     y = 1.0; */
/*   } */

/*   return y; */
/* } */

/*----------------*/
/*-Burn footprint-*/
/*----------------*/

/* Function: _trtt_3D_burn_footprint_3D
 *
 * Description:
 *
 *      This function calculates the impact of a given voxel and its footprint
 *      (integration) (B-spline model based) on the detector pixels, in the
 *      direct projection sense.
 *
 *      The footprint is assumed to be a 2D separable B-spline of degree 3,
 *      scaled by a magnification factor, and u-axis and v-axis distorsion
 *      factors (all embedded in a weighting factor normalization factors for
 *      coordinates on the detector plane).
 *
 *      >      |------|------|------|------|             ^
 *      >      |      |      |      |      |             |
 *      >      |      |      |      |      |             |
 *      >      |------|------|------|------|--- v_right  |
 *      >      |     _|______|___   |      |             |
 *      >      |    / |      |   \  |      |             |
 *      >      |---/--|------|----\-|------|             | nv
 *      >      |  |   |      |     ||      |             |
 *      >      |  |   |   o--|-----||------|--- vc       |
 *      >      |--|---|---'--|-----||------|             |
 *      >      |   \  |   '  |    / |      |             |
 *      >      |    \_|___'__|___/  |      |             |
 *      >      |------|---'--|------|------|--- v_left   v
 *      >      '          '         '
 *      >      '          '         '
 *      >      '          '         ' 
 *      >    u_left       uc      u_right
 *      >
 *      >      <--------------------------->
 *      >                   nu
 *
 *      The footprint being separable, the integration of the footprint is
 *      splitted into 2 mondimensional integrations. The mondimensional
 *      integration is calculated as the difference between the 2 cumulative
 *      integrations of the B-spline to the 2 detector pixel bounds normalized
 *      to the B-spline footprint frame (see
 *      <_trtt_cubic_bspline_integrale_3D>).
 *
 * Parameters: 
 *
 *      pix - data sinogram array to fill.
 *
 *      vox - value of the current B-spline coefficient (current voxel value)
 *
 *      uc - u-axis coordinate of the projection of the voxel center on the
 *           detector (homogeneous to detector pixels indexes).
 *
 *      vc - v-axis coordinate of the projection of the voxel center on the
 *           detector (homogeneous to detector pixels indexes).
 *
 *      u_left - u-axis index of the lower left corner detector pixel bound
 *               impinged by the footprint.
 *
 *      u_right - u-axis index of the upper right corner detector pixel bound
 *                impinged by the footprint.
 *
 *      v_left - v-axis index of the lower left corner detector pixel bound
 *               impinged by the footprint.
 *
 *      v_right - v-axis index of the upper right corner detector pixel bound
 *                impinged by the footprint.
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      rho0 - weighting factor for the footprint integration.
 *
 *      scl_norm_u - u-axis normalization factor for writting coordinates in the
 *                   centered normalized B-spline footprint frame.
 *
 *      scl_norm_v - v-axis normalization factor for writting coordinates in the
 *                   centered normalized B-spline footprint frame.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <_trtt_3D_burn_footprint_3D_transp>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
static int _trtt_3D_burn_footprint_3D(double *pix, 
				      const double vox,
				      const double uc,
				      const double vc,
				      const double u_left,
				      const double u_right,
				      const double v_left,
				      const double v_right,
				      const long nu,
				      const double rho0,
				      const double h_norm_u,
				      const double h_norm_v,
				      double scl_norm_u,
				      double scl_norm_v,
				      const double *footprint,
				      const double *cum_footprint,
				      const long size_footprint,
				      const double step_footprint,
				      const double *coord_footprint)
{
    long l, m;
    long lu_first, lu_next, mv_first, mv_next;
    double pu_first, u_next, pu_next;
    double pv_first, v_next, pv_next;    
    long npix_u, npix_v;
    double coord_footprint_start_u;
    double coord_footprint_start_v;
    double *c_integ_u;
    double c_integ_v;
    
    long bottom_left_corner_id = (long)v_left*nu + (long)u_left;
    long pix_id;
    // FIXME: we are in bound detector pixel coordinates
    // FIXME: 0 ; 1 ; 2 ; ... ; nu-1 ; nu
    // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv

    scl_norm_u = scl_norm_u/step_footprint;
    scl_norm_v = scl_norm_v/step_footprint;
    
    if (u_left<u_right && v_left<v_right) {

	/* Get the position of left-hand side of the footprint in pixel detector
	 * space */
	coord_footprint_start_u = ((u_left-uc)-h_norm_u*coord_footprint[0])*scl_norm_u;
	coord_footprint_start_v = ((v_left-vc)-h_norm_v*coord_footprint[0])*scl_norm_v;

	// FIXME: generic footprint coordinates for cubic b-spline
	// FIXME: -2 ; -1 ; 0 ; 1 ; 2
	npix_u = (long)(u_right - u_left);
	npix_v = (long)(v_right - v_left);

	/* Allocate "impact" coefficients in direction u */
	c_integ_u = TRTT_NEW_ARRAY0(double, npix_u);
	if (c_integ_u==NULL) return TRTT_FAILURE;
	
	/* Calculate "impact" coefficients in direction u */
	/* pu_first = _trtt_cubic_bspline_integrale_3D(coord_footprint_start_u); */
	lu_first = (long)floor(TRTT_MIN(TRTT_MAX(coord_footprint_start_u,0),size_footprint-1));
	pu_first = cum_footprint[lu_first];

	for (l=0; l<npix_u; ++l) { // for each impinged detector pixel in direction u
            /* Next detector pixel */
	    u_next = coord_footprint_start_u + (l+1)*scl_norm_u;
	    lu_next = (long)floor(TRTT_MIN(TRTT_MAX(u_next,0),size_footprint-1));
            /* Cumulative integration to the pixel bound */
	    /* pu_next = _trtt_cubic_bspline_integrale_3D(u_next); */
	    pu_next = cum_footprint[lu_next];
	    /* Calculate the impact */
	    c_integ_u[l] = pu_next-pu_first;
	    pu_first = pu_next;
	}
	
	/* Calculate "impact" coefficients in direction v */
	/* pv_first = _trtt_cubic_bspline_integrale_3D(coord_footprint_start_v); */
	mv_first = (long)floor(TRTT_MIN(TRTT_MAX(coord_footprint_start_v,0),size_footprint-1));
	pv_first = cum_footprint[mv_first];
	
	for (m=0; m<npix_v; ++m) { // for each impinged detector pixel in direction v
	    /* Next detector pixel */
	    v_next = coord_footprint_start_v + (m+1)*scl_norm_v;
	    mv_next = (long)floor(TRTT_MIN(TRTT_MAX(v_next,0),size_footprint-1));
	    /* Cumulative integration to the pixel bound */
	    /* pv_next = _trtt_cubic_bspline_integrale_3D(v_next); */
	    pv_next = cum_footprint[mv_next];
	    /* Find the "3D" index of the impinged detector pixel */
	    pix_id = bottom_left_corner_id + m*nu;
	    /* Calculate the impact */
	    c_integ_v = pv_next-pv_first;
	    for (l=0; l<npix_u; ++l) { // for each impinged detector pixel in direction u
		/* Accumulate in current detector pixel */
		*(pix + pix_id + l) += c_integ_v*c_integ_u[l]*vox*rho0;
	    }
	    pv_first = pv_next;
	}
	free(c_integ_u);
    }
    
    return TRTT_SUCCESS;
}

/* Function: _trtt_3D_DD_burn_footprint_3D
 *
 * Description:
 *
 *      This function calculates the impact of a given voxel and its footprint
 *      (integration) (Distance Driven model based) on the detector pixels in
 *      the direct projection sense. This function is only dedicated to the
 *      parallel beam geometry.
 *
 *      The Distance Driven model approximates the footprint of a cubic voxel by
 *      a rectangle on the detector plane. The footprint is obtained by
 *      projecting the corners of the central section of the cubic voxel.
 *
 *      >      |------|------|------|------|             ^
 *      >      |      |      |      |      |             |
 *      >      |      |      |      |      |             |
 *      >      |------|------|------|------|             |
 *      >      |     _|______|___ ..|......|....vb2      |
 *      >      |    | |      |   |  |      |             |
 *      >      |----|-|------|---|--|------|             | nv
 *      >      |    | |      |   |  |      |             |
 *      >      |    | |      |   |  |      |             |
 *      >      |----|-|------|---|--|------|             |
 *      >      |    | |      |   |  |      |             |
 *      >      |    |_|______|___|..|......|....vb1      |
 *      >      |------|------|------|------|             v
 *      >           '            '
 *      >           '            '
 *      >           '            ' 
 *      >          ub1          ub2
 *      >
 *      >      <--------------------------->
 *      >                   nu
 *
 *      The footprint being separable, the integration of the footprint is
 *      splitted into 2 mondimensional integrations.
 *
 * Parameters: 
 *
 *      pix - data sinogram array to fill.
 *
 *      vox - value of the current voxel.
 *
 *      ub1 - u-axis coordinate of the lower left corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      ub2 - u-axis coordinate of the upper right corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      vb1 - v-axis coordinate of the lower left corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      vb2 - v-axis coordinate of the upper right corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      u_scl - sampling step of detector pixels in u-direction.
 *
 *      v_scl - sampling step of detector pixels in v-direction.
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      h_norm - weighting factor for the footprint integration.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <_trtt_3D_burn_footprint_3D>, <_trtt_3D_DD_burn_footprint_3D_transp>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
static int _trtt_3D_DD_burn_footprint_3D(double *pix, 
					 const double vox,
					 const double ub1,
					 const double ub2,
					 const double vb1,
					 const double vb2,
					 const double u_scl,
					 const double v_scl,
					 const long nu,
					 const double h_norm)
{
    long l, m, i;
    double vup, vdown, uup, udown;
    double u_left, v_left;
    long npixu, npixu_eff, npixv, npixv_eff;
    double *c_integ_u, *c_integ_v;
    long bottom_left_corner_id, pix_id;

    /* Get the indexes of the lower left and upper right corners detector pixel
     * bounds impinged by the footprint */
    uup = ceil(ub1);             vup = ceil(vb1);
    udown = floor(ub2);          vdown = floor(vb2);
    /* Get the bounds positions of the footprint */
    u_left = floor(ub1);         v_left = floor(vb1);

    /* Find the "3D" index of the lower left impinged detector pixel */
    bottom_left_corner_id = (long)v_left*nu + (long)u_left;
    
    /* Find the number of impinged detector pixels */
    npixu = udown-uup;
    npixu_eff = ceil(ub2)-u_left;
    npixv = vdown-vup;
    npixv_eff = ceil(vb2)-v_left;

    if (npixu_eff>0.0 && npixv_eff>0.0) {
	
	/* Allocate "impact" coefficients in directions u and v */
	c_integ_u = TRTT_NEW_ARRAY0(double, npixu_eff);
	if (c_integ_u==NULL) return TRTT_FAILURE;
	c_integ_v = TRTT_NEW_ARRAY0(double, npixv_eff);
	if (c_integ_v==NULL) return TRTT_FAILURE;
	
	/* Calculate "impact" coefficients in direction u, "scanning" the
	 * impinged detecotr pixels and dealing with the different position
	 * cases */
	if (npixu>=0.0) {
	    l=0;
	    // for each impinged detector pixel in direction u
	    
	    // The pixel straddles the inferior bound
	    if(uup>ub1) {
		/* Calculate the impact */
		c_integ_u[l] = u_scl*(uup-ub1);
		l++;
	    }
	    // The pixel is entirely in the footprint
	    for (i=uup; i<udown; ++i) {
		/* Calculate the impact */
		c_integ_u[l] = u_scl;
		l++;
	    }
	    // The pixel straddles the superior bound
	    if (udown<ub2) {
		/* Calculate the impact */
		c_integ_u[l] = u_scl*(ub2-udown);
	    }
	} else {
	    // The lonely impinged pixel is entirely in the footprint

	    /* Calculate the impact */
	    c_integ_u[0] = u_scl*(ub2-ub1);
	}
	
	/* Calculate "impact" coefficients in direction v, "scanning" the
	 * impinged detecotr pixels and dealing with the different position
	 * cases */
	if (npixv>=0.0) {
	    m=0;
	    // for each impinged detector pixel in direction v

	    // The pixel straddles the inferior bound
	    if(vup>vb1) {
		/* Calculate the impact */
		c_integ_v[m] = v_scl*(vup-vb1);
		m++;
	    }

	    // The pixel is entirely in the footprint
	    for (i=vup; i<vdown; ++i) {
		/* Calculate the impact */
		c_integ_v[m] = v_scl;
		m++;
	    }
	    
	    // The pixel straddles the superior bound
	    if (vdown<vb2) {
		/* Calculate the impact */
		c_integ_v[m] = v_scl*(vb2-vdown);
	    }
	} else {
	    // The lonely impinged pixel is entirely in the footprint

	    /* Calculate the impact */
	    c_integ_v[0] = v_scl*(vb2-vb1);
	}
	
	
	for (m=0; m<npixv_eff; ++m) { // for each impinged detector pixel in direction v
	    pix_id = bottom_left_corner_id + m*nu;
	    for (l=0; l<npixu_eff; ++l) { // for each impinged detector pixel in direction u
		/* Accumulate in current detector pixel */
		*(pix + pix_id + l) += h_norm*c_integ_v[m]*c_integ_u[l]*vox;
	    }
	}

	free(c_integ_u);
	free(c_integ_v);
    }
    
    return TRTT_SUCCESS;
}
/*----------------------------*/
/*-Burn footprint (transpose)-*/
/*----------------------------*/

/* Function: _trtt_3D_burn_footprint_3D_transp
 *
 * Description:
 *
 *      This function calculates the impact of a given voxel and its footprint
 *      (integration) (B-spline model based) on the detector pixels, in the
 *      transpose projection sense.
 *
 *      The footprint is assumed to be a 2D separable B-spline of degree 3,
 *      scaled by a magnification factor, and u-axis and v-axis distorsion
 *      factors (all embedded in a weighting factor normalization factors for
 *      coordinates on the detector plane).
 *
 *      >      |------|------|------|------|             ^
 *      >      |      |      |      |      |             |
 *      >      |      |      |      |      |             |
 *      >      |------|------|------|------|--- v_right  |
 *      >      |     _|______|___   |      |             |
 *      >      |    / |      |   \  |      |             |
 *      >      |---/--|------|----\-|------|             | nv
 *      >      |  |   |      |     ||      |             |
 *      >      |  |   |   o--|-----||------|--- vc       |
 *      >      |--|---|---'--|-----||------|             |
 *      >      |   \  |   '  |    / |      |             |
 *      >      |    \_|___'__|___/  |      |             |
 *      >      |------|---'--|------|------|--- v_left   v
 *      >      '          '         '
 *      >      '          '         '
 *      >      '          '         ' 
 *      >    u_left       uc      u_right
 *      >
 *      >      <--------------------------->
 *      >                   nu
 *
 *      The footprint being separable, the integration of the footprint is
 *      splitted into 2 mondimensional integrations. The mondimensional
 *      integration is calculated as the difference between the 2 cumulative
 *      integrations of the B-spline to the 2 detector pixel bounds normalized
 *      to the B-spline footprint frame (see
 *      <_trtt_cubic_bspline_integrale_3D>).
 *
 * Parameters: 
 *
 *      vox - voxels image array.
 *
 *      pix - data sinogram array.
 *
 *      uc - u-axis coordinate of the projection of the voxel center on the
 *           detector (homogeneous to detector pixels indexes).
 *
 *      vc - v-axis coordinate of the projection of the voxel center on the
 *           detector (homogeneous to detector pixels indexes).
 *
 *      u_left - u-axis index of the lower left corner detector pixel bound
 *               impinged by the footprint.
 *
 *      u_right - u-axis index of the upper right corner detector pixel bound
 *                impinged by the footprint.
 *
 *      v_left - v-axis index of the lower left corner detector pixel bound
 *               impinged by the footprint.
 *
 *      v_right - v-axis index of the upper right corner detector pixel bound
 *                impinged by the footprint.
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      nx - size of the voxels image in x-direction.
 *
 *      ny - size of the voxels image in y-direction.
 *
 *      id_x - x-axis index of the current voxel to update.
 *
 *      id_y - y-axis index of the current voxel to update.
 *
 *      id_z - z-axis index of the current voxel to update.
 *
 *      rho0 - weighting factor for the footprint integration.
 *
 *      scl_norm_u - u-axis normalization factor for writting coordinates in the
 *                   centered normalized B-spline footprint frame.
 *
 *      scl_norm_v - v-axis normalization factor for writting coordinates in the
 *                   centered normalized B-spline footprint frame.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <_trtt_3D_burn_footprint_3D>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
static int _trtt_3D_burn_footprint_3D_transp(double *vox,
					     const double *pix,  
					     const double uc,
					     const double vc,
					     const double u_left,
					     const double u_right,
					     const double v_left,
					     const double v_right,
					     const long nu,
					     const long nx,
					     const long ny,
					     const long id_x,
					     const long id_y,
					     const long id_z,
					     const double rho0,
					     const double h_norm_u,
					     const double h_norm_v,
					     double scl_norm_u,
					     double scl_norm_v,
					     const double *footprint,
					     const double *cum_footprint,
					     const long size_footprint,
					     const double step_footprint,
					     const double *coord_footprint)
{
    long l, m;
    long lu_first, lu_next, mv_first, mv_next;
    double pu_first, u_next, pu_next;
    double pv_first, v_next, pv_next;    
    long npix_u, npix_v;
    double coord_footprint_start_u;
    double coord_footprint_start_v;
    double *c_integ_u;
    double c_integ_v;

    double *current_voxel = vox + (id_z*ny + id_y)*nx + id_x;
    
    long bottom_left_corner_id = (long)v_left * nu + (long)u_left;
    long pix_id;
    // FIXME: we are in bound detector pixel coordinates
    // FIXME: 0 ; 1 ; 2 ; ... ; nu-1 ; nu
    // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv

    scl_norm_u = scl_norm_u/step_footprint;
    scl_norm_v = scl_norm_v/step_footprint;
    
    if (u_left<u_right && v_left<v_right) {

	/* Get the position of left-hand side of the footprint in pixel detector
	 * space */
	coord_footprint_start_u = ((u_left-uc)-h_norm_u*coord_footprint[0])*scl_norm_u;
	coord_footprint_start_v = ((v_left-vc)-h_norm_v*coord_footprint[0])*scl_norm_v;
	// FIXME: generic footprint coordinates for cubic b-spline
	// FIXME: -2 ; -1 ; 0 ; 1 ; 2
	npix_u = (long)(u_right - u_left);
	npix_v = (long)(v_right - v_left);
	
	/* Allocate "impact" coefficients in direction u */
	c_integ_u = TRTT_NEW_ARRAY0(double, npix_u);
	if (c_integ_u==NULL) return TRTT_FAILURE;
	
	/* Calculate "impact" coefficients in direction u */
	/* pu_first = _trtt_cubic_bspline_integrale_3D(coord_footprint_start_u); */
	lu_first = (long)floor(TRTT_MIN(TRTT_MAX(coord_footprint_start_u,0),size_footprint-1));
	pu_first = cum_footprint[lu_first];
	
	for (l=0; l<npix_u; ++l) { // for each impinged detector pixel in direction u
	    /* Next detector pixel */
	    u_next = coord_footprint_start_u + (l+1)*scl_norm_u;
	    lu_next = (long)floor(TRTT_MIN(TRTT_MAX(u_next,0),size_footprint-1));
	    /* Cumulative integration to the pixel bound */
	    /* pu_next = _trtt_cubic_bspline_integrale_3D(u_next); */
	    pu_next = cum_footprint[lu_next];
	    /* Calculate the impact */
	    c_integ_u[l] = pu_next - pu_first;
	    pu_first = pu_next;
	}

	/* Calculate "impact" coefficients in direction v */
	/* pv_first = _trtt_cubic_bspline_integrale_3D(coord_footprint_start_v);	 */
	mv_first = (long)floor(TRTT_MIN(TRTT_MAX(coord_footprint_start_v,0),size_footprint-1));
	pv_first = cum_footprint[mv_first];
	
	for (m=0; m<npix_v; ++m) { // for each impinged detector pixel in direction v
	    v_next = coord_footprint_start_v + (m+1)*scl_norm_v;
	    mv_next = (long)floor(TRTT_MIN(TRTT_MAX(v_next,0),size_footprint-1));
	    /* Cumulative integration to the pixel bound */
	    /* pv_next = _trtt_cubic_bspline_integrale_3D(v_next); */
	    pv_next = cum_footprint[mv_next];
	    /* Find the "3D" index of the impinged detector pixel */
	    pix_id = bottom_left_corner_id + m*nu;
	    /* Calculate the impact */
	    c_integ_v = pv_next-pv_first;
	    for (l=0; l<npix_u; ++l) { // for each impinged detector pixel in direction u
		/* Accumulate in current detector pixel */
		*current_voxel += rho0 * c_integ_u[l] * c_integ_v * (*(pix + pix_id + l));
	    }
	    pv_first = pv_next;
	}
	free(c_integ_u);
    }
    return TRTT_SUCCESS;
}

/* Function: _trtt_3D_DD_burn_footprint_3D_transp
 *
 * Description:
 *
 *      This function calculates the impact of a given voxel and its footprint
 *      (integration) (Distance Driven model based) on the detector pixels in
 *      the transpose projection sense. This function is only dedicated to the
 *      parallel beam geometry.
 *
 *      The Distance Driven model approximates the footprint of a cubic voxel by
 *      a rectangle on the detector plane. The footprint is obtained by
 *      projecting the corners of the central section of the cubic voxel.
 *
 *      >      |------|------|------|------|             ^
 *      >      |      |      |      |      |             |
 *      >      |      |      |      |      |             |
 *      >      |------|------|------|------|             |
 *      >      |     _|______|___ ..|......|....vb2      |
 *      >      |    | |      |   |  |      |             |
 *      >      |----|-|------|---|--|------|             | nv
 *      >      |    | |      |   |  |      |             |
 *      >      |    | |      |   |  |      |             |
 *      >      |----|-|------|---|--|------|             |
 *      >      |    | |      |   |  |      |             |
 *      >      |    |_|______|___|..|......|....vb1      |
 *      >      |------|------|------|------|             v
 *      >           '            '
 *      >           '            '
 *      >           '            ' 
 *      >          ub1          ub2
 *      >
 *      >      <--------------------------->
 *      >                   nu
 *
 *      The footprint being separable, the integration of the footprint is
 *      splitted into 2 mondimensional integrations.
 *
 * Parameters: 
 *
 *      vox - voxels image array.
 *
 *      pix - data sinogram array.
 *
 *      ub1 - u-axis coordinate of the lower left corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      ub2 - u-axis coordinate of the upper right corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      vb1 - v-axis coordinate of the lower left corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      vb2 - v-axis coordinate of the upper right corner of the footprint
 *            (homogeneous to detector pixels indexes).
 *
 *      u_scl - sampling step of detector pixels in u-direction.
 *
 *      v_scl - sampling step of detector pixels in v-direction.
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      nx - size of the voxels image in x-direction.
 *
 *      ny - size of the voxels image in y-direction.
 *
 *      id_x - x-axis index of the current voxel to update.
 *
 *      id_y - y-axis index of the current voxel to update.
 *
 *      id_z - z-axis index of the current voxel to update.
 *
 *      h_norm - weighting factor for the footprint integration.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <_trtt_3D_burn_footprint_3D>, <_trtt_3D_DD_burn_footprint_3D_transp>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
static int _trtt_3D_DD_burn_footprint_3D_transp(double *vox, 
						const double *pix,
						const double ub1,
						const double ub2,
						const double vb1,
						const double vb2,
						const double u_scl,
						const double v_scl,
						const long nu,
						const long nx,
						const long ny,
						const long id_x,
						const long id_y,
						const long id_z,
						const double h_norm)
{
    long l, m, i;
    double vup, vdown, uup, udown;
    double u_left, v_left;
    long npixu, npixu_eff, npixv, npixv_eff;
    double *c_integ_u, *c_integ_v;
    long bottom_left_corner_id, pix_id;

    double *current_voxel = vox + (id_z*ny + id_y)*nx + id_x;

    /* Get the indexes of the lower left and upper right corners detector pixel
     * bounds impinged by the footprint */
    uup = ceil(ub1);             vup = ceil(vb1);
    udown = floor(ub2);          vdown = floor(vb2);
    /* Get the bounds positions of the footprint */
    u_left = floor(ub1);         v_left = floor(vb1);

    /* Find the "3D" index of the lower left impinged detector pixel */
    bottom_left_corner_id = (long)v_left*nu + (long)u_left;
    
    /* Find the number of impinged detector pixels */
    npixu = udown-uup;
    npixu_eff = ceil(ub2)-u_left;
    npixv = vdown-vup;
    npixv_eff = ceil(vb2)-v_left;

    if (npixu_eff>0.0 && npixv_eff>0.0) {
	
	/* Allocate "impact" coefficients in directions u and v */
	c_integ_u = TRTT_NEW_ARRAY0(double, npixu_eff);
	if (c_integ_u==NULL) return TRTT_FAILURE;
	c_integ_v = TRTT_NEW_ARRAY0(double, npixv_eff);
	if (c_integ_v==NULL) return TRTT_FAILURE;
	
	/* Calculate "impact" coefficients in direction u, "scanning" the
	 * impinged detecotr pixels and dealing with the different position
	 * cases */
	if (npixu>=0.0) {
	    l=0;
	    // for each impinged detector pixel in direction u

	    // The pixel straddles the inferior bound
	    if(uup>ub1) {
		/* Calculate the impact */
		c_integ_u[l] = u_scl*(uup-ub1);
		l++;
	    }
	    // The pixel is entirely in the footprint
	    for (i=uup; i<udown; ++i) {
		/* Calculate the impact */
		c_integ_u[l] = u_scl;
		l++;
	    }
	    // The pixel straddles the superior bound
	    if (udown<ub2) {
		/* Calculate the impact */
		c_integ_u[l] = u_scl*(ub2-udown);
	    }
	} else {
	    // The lonely impinged pixel is entirely in the footprint
	    
	    /* Calculate the impact */	    
	    c_integ_u[0] = u_scl*(ub2-ub1);
	}
	
	/* Calculate "impact" coefficients in direction v, "scanning" the
	 * impinged detecotr pixels and dealing with the different position
	 * cases */
	if (npixv>=0.0) {
	    m=0;
	    // for each impinged detector pixel in direction v

	    // The pixel straddles the inferior bound
	    if(vup>vb1) {
		/* Calculate the impact */
		c_integ_v[m] = v_scl*(vup-vb1);
		m++;
	    }

	    // The pixel is entirely in the footprint
	    for (i=vup; i<vdown; ++i) {
		/* Calculate the impact */
		c_integ_v[m] = v_scl;
		l++;
	    }

	    // The pixel straddles the superior bound
	    if (vdown<vb2) {
		/* Calculate the impact */
		c_integ_v[m] = v_scl*(vb2-vdown);
	    }
	} else {
	    // The lonely impinged pixel is entirely in the footprint

	    /* Calculate the impact */
	    c_integ_v[0] = v_scl*(vb2-vb1);
	}
	
	
	for (m=0; m<npixv_eff; ++m) { // for each impinged detector pixel in direction v
	    pix_id = bottom_left_corner_id + m*nu;
	    for (l=0; l<npixu_eff; ++l) { // for each impinged detector pixel in direction u
		/* Accumulate in current detector pixel */
		*current_voxel += h_norm*c_integ_v[m]*c_integ_u[l]*(*(pix + pix_id + l));
	    }
	}

	free(c_integ_u);
	free(c_integ_v);
    }
    
    return TRTT_SUCCESS;
}

/* PUBLIC FUNCTIONS =================================================== */

/*---------------------------------------------------------------*/
/* See NaturalDocs Documentation for the help of these functions */
/*---------------------------------------------------------------*/

//////////////////
// 3D OPERATORS //
//////////////////

/*------------------*/
/*-3D PARALLEL BEAM-*/
/*------------------*/
int trtt_PB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job)
{
    /*  Calculated parameters */
    double in_ijk;
    double uc, u_left, u_right, vc, v_left, v_right;
    double uc_temp1, vc_temp1;
    double uc_temp2, vc_temp2;

    long i, j, k;
    
    /* Extract some pre-calculated coefficients */
    double a1=xsd[0]; double a2=xsd[1]; double a3=xsd[2]; double a4=xsd[3];
    double a5=xsd[4]; double a6=xsd[5]; double a7=xsd[6]; double a8=xsd[7];
    double a9=xsd[8]; double a10=xsd[9]; double a11=xsd[10]; double a12=xsd[11];
    
    /* double o1=xos[0]; double o2=xos[1]; double o3=xos[2]; double o4=xos[3]; */
    double o5=xos[4]; double o6=xos[5]; double o7=xos[6]; double o8=xos[7];
    double o9=xos[8]; double o10=xos[9]; double o11=xos[10]; double o12=xos[11];
    
    double a1inv = -1.0*(1.0/a1);
    double b1 = a6 + a1inv*a5*a2;      double c1 = a10 + a1inv*a9*a2;
    double b2 = a7 + a1inv*a5*a3;      double c2 = a11 + a1inv*a9*a3;
    double b3 = a8 + a1inv*a5*a4;      double c3 = a12 + a1inv*a9*a4;

    double v_x = b1*o5 + b2*o9;        double u_x = c1*o5 + c2*o9;
    double v_y = b1*o6 + b2*o10;       double u_y = c1*o6 + c2*o10;
    double v_z = b1*o7 + b2*o11;       double u_z = c1*o7 + c2*o11;
    double v_c = b1*o8 + b2*o12 + b3;  double u_c = c1*o8 + c2*o12 + c3;

    double a1sqr = a1*a1;
    double a5sqr = v_scl*v_scl*a5*a5; //FIXME: metric scale
    double a9sqr = u_scl*u_scl*a9*a9; //FIXME: metric scale

    /* Initialize scaling factors */
    double delta_u = sqrt(1.0+(a9sqr/(a1sqr+a5sqr)));
    double delta_v = sqrt(1.0+(a5sqr/a1sqr));
    double rho0 = s_scl*s_scl*s_scl*delta_u*delta_v;
    double h_norm_u = s_scl/u_scl;
    double h_norm_v = s_scl/v_scl;
    double scl_norm_u = 1.0/h_norm_u;
    double scl_norm_v = 1.0/h_norm_v;
    /* Retrieve the half support of the footprint on detector */
    double half_support_init = 0.5*(s_deg+1);
    double half_support_u = half_support_init*h_norm_u;
    double half_support_v = half_support_init*h_norm_v;
    
    for (k=0; k<nz; ++k) {

	/* Projection of the voxel center */
	vc_temp1 = v_z*k + v_c;
	uc_temp1 = u_z*k + u_c;

	for (j=0; j<ny; ++j) {

	    /* Projection of the voxel center */
	    vc_temp2 = vc_temp1 + v_y*j;
	    uc_temp2 = uc_temp1 + u_y*j;

	    for (i=0; i<nx; ++i) {
		
		/* Projection of the voxel center */
		vc = vc_temp2 + v_x*i;
		uc = uc_temp2 + u_x*i;

		/* Get the bounds positions of the footprint */
		u_left = uc-half_support_u;
		u_right = uc+half_support_u;
		v_left = vc-half_support_v;
		v_right = vc+half_support_v;
		
		/* Constrain the bounds positions to the detector support */
		if (u_left < 0.0) u_left = 0.0;
		else if (u_left >= 0.0 && u_left <= (nu-1)) u_left = floor(u_left);
		else u_left = nu-1;
            
		if (u_right < 0.0) u_right = 0.0;
		else if (u_right >= 0.0 && u_right <= (nu-1)) u_right = floor(u_right)+1;
		else u_right = nu;

		if (v_left < 0.0) v_left = 0.0;
		else if (v_left >= 0.0 && v_left <= (nv-1)) v_left = floor(v_left);
		else v_left = nv-1;
            
		if (v_right < 0.0) v_right = 0.0;
		else if (v_right >= 0.0 && v_right <= (nv-1)) v_right = floor(v_right)+1;
		else v_right = nv;

		/* Burn footprint on Detector */
		if (job == 0) {
		    /* DIRECT OPERATOR */
		    /* Get voxel value */
		    in_ijk = TRTT_VOXEL_VALUE_ACCESS(in,nx,ny,i,j,k);
		    if (_trtt_3D_burn_footprint_3D(out, in_ijk, uc, vc, u_left, u_right, \
						   v_left, v_right, nu, rho0, h_norm_u, h_norm_v, \
						   scl_norm_u, scl_norm_v, \
						   footprint, cum_footprint, size_footprint, \
						   step_footprint, coord_footprint) == TRTT_FAILURE) {
			return TRTT_FAILURE;
		    }  
		} else {
		    /* TRANSPOSE OPERATOR */
		    if (_trtt_3D_burn_footprint_3D_transp(out, in, uc, vc, u_left, u_right, \
							  v_left, v_right, nu, \
							  nx, ny, i, j, k, \
							  rho0, h_norm_u, h_norm_v, \
							  scl_norm_u, scl_norm_v, \
							  footprint, cum_footprint, size_footprint, \
							  step_footprint, coord_footprint) == TRTT_FAILURE) {
			return TRTT_FAILURE;
		    }
		}
	    }
	}	
    }

    return TRTT_SUCCESS;
}

int trtt_DD_PB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job)
{
    /*  Calculated parameters */
    double in_ijk, uc, vc;
    double uc_temp1, vc_temp1;
    double uc_temp2, vc_temp2;
    double vb1, vb2, ub1, ub2;
    
    long i, j, k;
    
    /* Extract some pre-calculated coefficients */
    double a1=xsd[0]; double a2=xsd[1]; double a3=xsd[2]; double a4=xsd[3];
    double a5=xsd[4]; double a6=xsd[5]; double a7=xsd[6]; double a8=xsd[7];
    double a9=xsd[8]; double a10=xsd[9]; double a11=xsd[10]; double a12=xsd[11];
    
    /* double o1=xos[0]; double o2=xos[1]; double o3=xos[2]; double o4=xos[3]; */
    double o5=xos[4]; double o6=xos[5]; double o7=xos[6]; double o8=xos[7];
    double o9=xos[8]; double o10=xos[9]; double o11=xos[10]; double o12=xos[11];
    
    double a1inv = -1.0*(1.0/a1);
    double b1 = a6 + a1inv*a5*a2;      double c1 = a10 + a1inv*a9*a2;
    double b2 = a7 + a1inv*a5*a3;      double c2 = a11 + a1inv*a9*a3;
    double b3 = a8 + a1inv*a5*a4;      double c3 = a12 + a1inv*a9*a4;

    double v_x = b1*o5 + b2*o9;        double u_x = c1*o5 + c2*o9;
    double v_y = b1*o6 + b2*o10;       double u_y = c1*o6 + c2*o10;
    double v_z = b1*o7 + b2*o11;       double u_z = c1*o7 + c2*o11;
    double v_c = b1*o8 + b2*o12 + b3;  double u_c = c1*o8 + c2*o12 + c3;

    double Du = TRTT_ABS(u_z);
    double Dv;

    double h_norm;

    /* Determine voxel section direction (x or y) */
    if (TRTT_ABS(v_y)>TRTT_ABS(v_x)) Dv = TRTT_ABS(v_y);
    else                             Dv = TRTT_ABS(v_x);

    h_norm = (s_scl*s_scl*s_scl)/(u_scl*v_scl*Du*Dv);

    Du = 0.5*Du;
    Dv = 0.5*Dv;
    
    for (k=0; k<nz; ++k) {
	/* Projection of the voxel center */
	vc_temp1 = v_z*k + v_c;
	uc_temp1 = u_z*k + u_c;
	
	for (j=0; j<ny; ++j) {
	    /* Projection of the voxel center */
	    vc_temp2 = vc_temp1 + v_y*j;
	    uc_temp2 = uc_temp1 + u_y*j;
	    
	    for (i=0; i<nx; ++i) {
		/* Projection of the voxel center */
		vc = vc_temp2 + v_x*i;
		uc = uc_temp2 + u_x*i;

		/* Determine the corners positions of the rectangular
		 * footprint */
		ub1 = uc-Du;                 vb1 = vc-Dv;
		ub2 = uc+Du;                 vb2 = vc+Dv;
		
                /* Constrain the corners to the detector support */
		ub1 = TRTT_MIN(TRTT_MAX(ub1,0),nu-1);
		ub2 = TRTT_MIN(TRTT_MAX(ub2,0),nu-1);

		vb1 = TRTT_MIN(TRTT_MAX(vb1,0),nv-1);
		vb2 = TRTT_MIN(TRTT_MAX(vb2,0),nv-1);

                /* Burn footprint on Detector */
		if (job == 0) {
		    /* DIRECT OPERATOR */
		    /* Get voxel value */
		    in_ijk = TRTT_VOXEL_VALUE_ACCESS(in,nx,ny,i,j,k);
		    if (_trtt_3D_DD_burn_footprint_3D(out, in_ijk, \
						      ub1, ub2, vb1, vb2, \
						      u_scl, v_scl, nu, \
						      h_norm) == TRTT_FAILURE) {
			return TRTT_FAILURE;
		    }
		} else {
		    /* TRANSPOSE OPERATOR */
		    if (_trtt_3D_DD_burn_footprint_3D_transp(out, in,	\
							     ub1, ub2, vb1, vb2, \
							     u_scl, v_scl, nu, \
							     nx, ny, i, j, k, \
							     h_norm) == TRTT_FAILURE) {
			return TRTT_FAILURE;
		    }
		}
	    }
	}
    }
    
    return TRTT_SUCCESS;
}


/*--------------*/
/*-3D CONE BEAM-*/
/*--------------*/
int trtt_CB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job)
{
    /*  Calculated parameters */
    double in_ijk;
    double num_u, num_v, den, wk;
    double num_u_temp1, num_v_temp1, den_temp1, wk_temp1;
    double num_u_temp2, num_v_temp2, den_temp2, wk_temp2;
    double uc, u_left, u_right, vc, v_left, v_right;

    long i, j, k;
    
    /* Extract some pre-calculated coefficients */
    double a1=xsd[0]; double a2=xsd[1]; double a3=xsd[2]; double a4=xsd[3];
    double a5=xsd[4]; double a6=xsd[5]; double a7=xsd[6]; double a8=xsd[7];
    double a9=xsd[8]; double a10=xsd[9]; double a11=xsd[10]; double a12=xsd[11];
    
    double o1=xos[0]; double o2=xos[1]; double o3=xos[2]; double o4=xos[3];
    double o5=xos[4]; double o6=xos[5]; double o7=xos[6]; double o8=xos[7];
    double o9=xos[8]; double o10=xos[9]; double o11=xos[10]; double o12=xos[11];
    
    double a1inv = -1.0*(1.0/a1);
    double b1 = a6 + a1inv*a5*a2; double c1 = a10 + a1inv*a9*a2;
    double b2 = a7 + a1inv*a5*a3; double c2 = a11 + a1inv*a9*a3;
    double b3 = a8 + a1inv*a5*a4; double c3 = a12 + a1inv*a9*a4;

    double num_v_x = b1*o5 + b2*o9;   double num_u_x = c1*o5 + c2*o9;
    double num_v_y = b1*o6 + b2*o10;  double num_u_y = c1*o6 + c2*o10;
    double num_v_z = b1*o7 + b2*o11;  double num_u_z = c1*o7 + c2*o11;
    double num_v_c = b1*o8 + b2*o12;  double num_u_c = c1*o8 + c2*o12;

    double den_x = a1*o1 + a2*o5 + a3*o9;
    double den_y = a1*o2 + a2*o6 + a3*o10;
    double den_z = a1*o3 + a2*o7 + a3*o11;
    double den_c = a1*o4 + a2*o8 + a3*o12;

    double a4sqr = a4*a4;
    double a4m = -1.0*a4;

    /* Various intermediate variables */
    double a4quot;
    double uca12, vca8, vca8sqr;

    /* Initialize scaling factors */
    double rho0;
    double rho0_init = s_scl*s_scl*s_scl;
    double delta_u, delta_v, magn, magn_sqr;
    /* Normalized b-spline step to pixel detector sampling */
    double h_norm_u, h_norm_v;
    double h_norm_u_init = s_scl/u_scl;
    double h_norm_v_init = s_scl/v_scl;
    double scl_norm_u, scl_norm_v;
    /* Retrieve the half support of the footprint on detector */
    double half_support_u, half_support_v;
    double half_support_init = 0.5*(s_deg+1);
    
    for (k=0; k<nz; ++k) {

	/* Projection of the voxel center */
	wk_temp1 = o3*k + o4;
	num_u_temp1 = num_u_z*k + num_u_c;
	num_v_temp1 = num_v_z*k + num_v_c;
	den_temp1 = den_z*k + den_c;

	for (j=0; j<ny; ++j) {

	    /* Projection of the voxel center */
	    wk_temp2 = wk_temp1 + o2*j;
	    num_u_temp2 = num_u_temp1 + num_u_y*j;
	    num_v_temp2 = num_v_temp1 + num_v_y*j;
	    den_temp2 = den_temp1 + den_y*j;

	    for (i=0; i<nx; ++i) {
		
		/* Projection of the voxel center */
		wk = wk_temp2 + o1*i;
		num_u = num_u_temp2 + num_u_x*i;
		num_v = num_v_temp2 + num_v_x*i;
		den = den_temp2 + den_x*i;

		/* Projection on Detector */
		a4quot = a4m/den;
		uc = a4quot*num_u + c3;
		vc = a4quot*num_v + b3;

		/* Calculate the scaling factors */
		if (wk != 0) {   
		    uca12 = (uc-a12)*u_scl; //FIXME: metric scale
		    vca8 = (vc-a8)*v_scl;   //FIXME: metric scale
		    vca8sqr = vca8*vca8;
		    magn = a4m/wk;
		    magn_sqr = magn*magn;
		    delta_u = sqrt(1.0+((uca12*uca12)/(a4sqr+vca8sqr)));
		    delta_v = sqrt(1.0+((vca8sqr)/(a4sqr)));
		    rho0 = rho0_init*magn_sqr*delta_u*delta_v;

		    h_norm_u = h_norm_u_init*magn*delta_u;
		    h_norm_v = h_norm_v_init*magn*delta_v;
		    scl_norm_u = 1.0/h_norm_u;
		    scl_norm_v = 1.0/h_norm_v;

		    half_support_u = half_support_init * h_norm_u;
		    half_support_v = half_support_init * h_norm_v;    
		} else {
		    fprintf(stderr, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n");
		    return TRTT_FAILURE;
		}

		/* Get the bounds positions of the footprint */
		u_left = uc-half_support_u;
		u_right = uc+half_support_u;
		v_left = vc-half_support_v;
		v_right = vc+half_support_v;

		/* Constrain the bounds positions to the detector support */
		if (u_left < 0.0) u_left = 0.0;
		else if (u_left >= 0.0 && u_left <= (nu-1)) u_left = floor(u_left);
		else u_left = nu-1;
            
		if (u_right < 0.0) u_right = 0.0;
		else if (u_right >= 0.0 && u_right <= (nu-1)) u_right = floor(u_right)+1;
		else u_right = nu;

		if (v_left < 0.0) v_left = 0.0;
		else if (v_left >= 0.0 && v_left <= (nv-1)) v_left = floor(v_left);
		else v_left = nv-1;
            
		if (v_right < 0.0) v_right = 0.0;
		else if (v_right >= 0.0 && v_right <= (nv-1)) v_right = floor(v_right)+1;
		else v_right = nv;

		/* Burn footprint on Detector */
		if (job == 0) {
		    /* DIRECT OPERATOR */
		    /* Get voxel value */
		    in_ijk = TRTT_VOXEL_VALUE_ACCESS(in,nx,ny,i,j,k);

		    if (_trtt_3D_burn_footprint_3D(out, in_ijk, uc, vc, u_left, u_right, \
						   v_left, v_right, nu, rho0, h_norm_u, h_norm_v, \
						   scl_norm_u, scl_norm_v, \
						   footprint, cum_footprint, size_footprint, \
						   step_footprint, coord_footprint) == TRTT_FAILURE) {
			return TRTT_FAILURE;
		    }  
		} else {
		    /* TRANSPOSE OPERATOR */
		    if (_trtt_3D_burn_footprint_3D_transp(out, in, uc, vc, u_left, u_right, \
							  v_left, v_right, nu, nx, ny, i, j, k, \
							  rho0, h_norm_u, h_norm_v, \
							  scl_norm_u, scl_norm_v, \
							  footprint, cum_footprint, size_footprint, \
							  step_footprint, coord_footprint) == TRTT_FAILURE) {
			return TRTT_FAILURE;
		    }
		}
	    }
	}	
    }

    return TRTT_SUCCESS;
}

/* int trtt_DD_CB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job) */
/* { */
/*     /\*  Calculated parameters *\/ */
/*     double in_ijk; */
/*     double num_u_bb, num_u_ub, num_u_bu, num_u_uu; */
/*     double num_v_bb, num_v_ub, num_v_bu, num_v_uu; */
/*     double den_bb, den_ub, den_bu, den_uu; */
/*     double wk, wk_temp1, wk_temp2; */
/*     double vb1, vb2, ub1, ub2; */
    
/*     long i, j, k; */
/*     int axis; */
    
/*     /\* Extract some pre-calculated coefficients *\/ */
/*     double a1=xsd[0]; double a2=xsd[1]; double a3=xsd[2]; double a4=xsd[3]; */
/*     double a5=xsd[4]; double a6=xsd[5]; double a7=xsd[6]; double a8=xsd[7]; */
/*     double a9=xsd[8]; double a10=xsd[9]; double a11=xsd[10]; double a12=xsd[11]; */
    
/*     double o1=xos[0]; double o2=xos[1]; double o3=xos[2]; double o4=xos[3]; */
/*     double o5=xos[4]; double o6=xos[5]; double o7=xos[6]; double o8=xos[7]; */
/*     double o9=xos[8]; double o10=xos[9]; double o11=xos[10]; double o12=xos[11]; */
    
/*     double a1inv = -1.0*(1.0/a1); */
/*     double b1 = a6 + a1inv*a5*a2;      double c1 = a10 + a1inv*a9*a2; */
/*     double b2 = a7 + a1inv*a5*a3;      double c2 = a11 + a1inv*a9*a3; */
/*     double b3 = a8 + a1inv*a5*a4;      double c3 = a12 + a1inv*a9*a4; */

/*     double num_v_x = b1*o5 + b2*o9;   double num_u_x = c1*o5 + c2*o9; */
/*     double num_v_y = b1*o6 + b2*o10;  double num_u_y = c1*o6 + c2*o10; */
/*     double num_v_z = b1*o7 + b2*o11;  double num_u_z = c1*o7 + c2*o11; */
/*     double num_v_c = b1*o8 + b2*o12;  double num_u_c = c1*o8 + c2*o12; */

/*     double den_x = a1*o1 + a2*o5 + a3*o9; */
/*     double den_y = a1*o2 + a2*o6 + a3*o10; */
/*     double den_z = a1*o3 + a2*o7 + a3*o11; */
/*     double den_c = a1*o4 + a2*o8 + a3*o12; */

/*     double a4sqr = a4*a4; */
/*     double a4m = -1.0*a4; */

/*     double h_norm; */
/*     double h_norm_init = s_scl*s_scl*s_scl; */

/*     /\* Determine voxel section direction (x or y) *\/ */
/*     if (TRTT_ABS(v_y)>TRTT_ABS(v_x)) axis = 2; */
/*     else                             axis = 1; */
    
/*     for (k=0; k<nz; ++k) { */


	
/* 	for (j=0; j<ny; ++j) { */


	    
/* 	    for (i=0; i<nx; ++i) { */



/* 		ub1 = uc-Du;                 vb1 = vc-Dv; */
/* 		ub2 = uc+Du;                 vb2 = vc+Dv; */

/* 		ub1 = TRTT_MIN(TRTT_MAX(ub1,0),nu-1); */
/* 		ub2 = TRTT_MIN(TRTT_MAX(ub2,0),nu-1); */

/* 		vb1 = TRTT_MIN(TRTT_MAX(vb1,0),nv-1); */
/* 		vb2 = TRTT_MIN(TRTT_MAX(vb2,0),nv-1); */

/*                 /\* Burn footprint on Detector *\/ */
/* 		if (job == 0) { */
/* 		    /\* DIRECT OPERATOR *\/ */
/* 		    /\* Get voxel value *\/ */
/* 		    in_ijk = TRTT_VOXEL_VALUE_ACCESS(in,nx,ny,i,j,k); */
/* 		    if (_trtt_3D_DD_burn_footprint_3D(out, in_ijk, \ */
/* 						      ub1, ub2, vb1, vb2, \ */
/* 						      u_scl, v_scl, nu, \ */
/* 						      h_norm) == TRTT_FAILURE) { */
/* 			return TRTT_FAILURE; */
/* 		    } */
/* 		} else { */
/* 		    /\* TRANSPOSE OPERATOR *\/ */
/* 		    if (_trtt_3D_DD_burn_footprint_3D_transp(out, in,	\ */
/* 							     ub1, ub2, vb1, vb2, \ */
/* 							     u_scl, v_scl, nu, \ */
/* 							     nx, ny, i, j, k, \ */
/* 							     h_norm) == TRTT_FAILURE) { */
/* 			return TRTT_FAILURE; */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */
/*     } */
    
/*     return TRTT_SUCCESS; */
/* } */
