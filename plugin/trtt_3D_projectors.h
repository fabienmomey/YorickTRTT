/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * trtt_3D_projectors.h --
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

/* PREPROCESSOR =============================================================== */

#ifndef TRTT_3D_PROJECTORS_H
#define TRTT_3D_PROJECTORS_H 1

/* PREPROCESSOR INCLUSIONS ==================================================== */

#include <stdio.h>
#include "trtt_common.h"
#include "trtt_xform3d.h"

/* PROTOYPES ================================================================== */

TRTT_BEGIN_C_DECLS

/*
 * Title: TRTT Radon projectors for tomography
 *
 *   Version 1.0
 *
 *   Copyright (C) 2011, the MiTiV project.
 *
 *   This package contains the projectors performing the algebraic Radon
 *   transform at several dimensions and geometries.
 *
 *   The definitions are available by: 
 *   > #include <trtt_3D_projectors.h>
 */

/*---------------------------------------------------------------------------*/

/*
 * Function: trtt_PB3D
 *
 *      3D projector based on the B-spline model, in parallel beam geometry.
 *
 * Description:
 *
 *      This function is the projector performing the 3D algebraic Radon
 *      transform in parallel beam geometry, for a given source/detector
 *      position (i.e. a single projection), using a data modelization based
 *      on a representation of the voxels image using B-splines instead of the
 *      classical "staircase" functions. 
 *
 *      The calibration parameters are passed in arguments of the function, so
 *      are the voxels image and the data pixels. The calibration parameters
 *      represents a single projection. In addition to parameters giving the
 *      dimensions and scales of the image and the data (see section
 *      Parameters), the projection geometry features are stored as 2
 *      transform matrices:
 *
 *      OBJECT TO SOURCE - *XOS* (for Object to Source) => the grid of voxels
 *                                                         coordinates are
 *                                                         rotated (*Ros*
 *                                                         transform below) in
 *                                                         a frame
 *                                                         corresponding to
 *                                                         the detector frame,
 *                                                         orthogonal to the
 *                                                         direction of rays
 *                                                         propagation, and
 *                                                         translated (*Tos*
 *                                                         transform below)
 *                                                         such that the
 *                                                         origin of the frame
 *                                                         is the source.
 *
 *      > XOS = | o1  o2  o3  o4  | = Tos . Ros . Xoi
 *      >       | o5  o6  o7  o8  |
 *      >       | o9  o10 o11 o12 |
 *      >       | 0   0   0   1   |
 *
 *      *XOS* also makes pass from the voxels indexes to the metric scale, the
 *      origin of which is taken at the center of the object image: in the
 *      projector, a loop on the voxels indexes (from i = 0 -> n-1 in each
 *      direction) is made => the conversion to the metric scale is therefore
 *      *embedded* in *XOS* (*Xoi* transform).
 *      
 *      > Xoi = | 1  0  0  -(1/2).(nx-1).s_scl |
 *      >       | 0  1  0  -(1/2).(ny-1).s_scl |
 *      >       | 0  0  1  -(1/2).(nz-1).s_scl |
 *      >       | 0  0  0             1        |
 *      >
 *      > Centers of the voxels are taken as references for indexes. 
 *      > They are indexed from 0 to nx-1.
 *
 *      > Tos = | 1  0  0  Rsc |
 *      >       | 0  1  0   0  |     where Rsc is the source
 *      >       | 0  0  1   0  |     to object center distance.
 *      >       | 0  0  0   1  |
 *
 *
 *      *Ros* is specified by a *quaternion* *qw*, *qx*, *qy* and *qz*. A
 *      rotation of angle *theta* counterclockwise around an axis of normalized
 *      direction vector (*ex*, *ey*, *ez*) is defined by this quaternion.
 *
 *      o *qw* = cos(theta/2)       ;   QW = qw.qw
 *      o *qx* = ex.sin(theta/2)  ;   QX = qx.qx
 *      o *qy* = ey.sin(theta/2)  ;   QY = qy.qy
 *      o *qz* = ez.sin(theta/2)  ;   QZ = qz.qz
 *
 *      > Ros = | QW-QX-QY-QZ      2.qx.qy-2.qw.qz  2.qx.qz+2.qw.qy  0 |
 *      >       | 2.qx.qy+2.qw.qz  QW-QX+QY-QZ      2.qy.qz-2.qw.qx  0 |
 *      >       | 2.qx.qz-2.qw.qy  2.qy.qz+2.qw.qx  QW-QX-QY+QZ      0 |
 *      >       | 0                0                0                1 |
 *
 *      Example - Case of a *C-arm scanner*, where the source
 *      simply rotates around the z-axis
 *
 *      | o ex = ey = 0
 *      | o ez = 1
 *      |
 *      | Ros = | cos(theta)      -sin(theta)  0  0 |
 *      |       | sin(theta)      cos(theta)   0  0 |
 *      |       | 0               0            0  0 |
 *      |       | 0               0            0  1 |
 *      |
 *      | where theta is the angle from the w-axis, orthogonal to the detector
 *      | plane, and the x-axis of the object frame.
 *      |
 *      | It is no more than a 3D rotation of theta angle around one of the
 *      | 3 cartesian axis.
 *
 *      SOURCE TO DETECTOR - *XSD* (for Source to Detector) => the
 *                                                             coefficients of
 *                                                             this transform
 *                                                             are used to
 *                                                             project the
 *                                                             voxels
 *                                                             positions on
 *                                                             the detector
 *                                                             plane,
 *                                                             according to
 *                                                             geometry of
 *                                                             projection. This
 *                                                             matrix contains
 *                                                             for instance
 *                                                             the focale
 *                                                             distance from
 *                                                             the source to
 *                                                             the detector
 *                                                             plane, as well
 *                                                             as the offsets
 *                                                             between the
 *                                                             center of the
 *                                                             detector and
 *                                                             the orthogonal
 *                                                             projection of
 *                                                             the source
 *                                                             point on the
 *                                                             detector plane.
 *
 *      > XSD = | a1  a2  a3  a4  | = Xdi^(-1) . Tsd . Rsd
 *      >       | a5  a6  a7  a8  |
 *      >       | a9  a10 a11 a12 |
 *      >       | 0   0   0   1   |
 *
 *      *XSD* also makes pass from the metric scale of the detector, the
 *      origin of which is taken at the center of the detector plane, to
 *      detector pixels indexes (from i = 0 -> n in each direction) => the
 *      detector pixels scale normalization is therefore *embedded* in *XSD*
 *      (*Xsi* transform).
 *      
 *      > Xsi = | 1  0  0         1        |
 *      >       | 0  1  0  -(1/2).nv.v_scl |
 *      >       | 0  0  1  -(1/2).nu.u_scl |
 *      >       | 0  0  0         1        |
 *      >
 *      > Detector pixels bounds are taken as references for indexes. 
 *      > They are indexed from 0 to nu [nv respectively].
 *
 *      > Tsd = | 1  0  0     -f    |     where f is the focale
 *      >       | 0  1  0   -v_off  |           v_off is the offset in on v-axis
 *      >       | 0  0  1   -u_off  |           u_off is the offset in on u-axis
 *      >       | 0  0  0      1    |
 *      >
 *      > The offsets are the distances on each axis from the orthogonal
 *      > projection of the source and the center of the detector plane.
 *
 *      The same way as *Ros*, *Rsd* is specified by a *quaternion* *qw*,
 *      *qx*, *qy* and *qz*.
 *
 *      Important notes - Distinguish geometries:
 *      o PARALLEL BEAM: In most standard cases, *Rsd* will be equal to the
 *      identity matrix, because the detector is almost always assumed to be
 *      orthogonal to the direction of rays propagation. However it is
 *      possible to suppose a tip/tilt of the detector plane from this typicla
 *      direction, which can be implemented in this transform *Rsd*.
 *      o CONE BEAM: The previous case does not hold for this geometry,
 *      because we have a divergent beam. As a result a tip/tilt of the
 *      detector plane can be directly implemented in the first frame rotation
 *      *Ros*.
 *
 *      REMARK - *XOS* and *XSD* are 4x4 matrices to be in homogeneous
 *      coordinate system to allow translations. Then in practice only 12
 *      coefficients are necessary and are extracted in the projector
 *      implementation. As a result the parameters <xos> and <xsd> are simply
 *      arrays of 12 elements, numbered as shown above.
 *
 *      HOW TO CREATE THE XOS AND XSD COEFFICIENTS - Use the tools of XFORM3D
 *      (definitions in <trtt_xform3d.h>) to create 3D coordinate transform.
 *      o Quaternion rotations.
 *      o Translations.
 *      o Scaling.
 *      o Combination, inversion of matrices.
 *      o Retrieve coefficients.   
 *
 * Parameters:
 *
 *      out - output array (must be declared first, with right dimensions):
 *            it can be either the data sinogram if the direct projector is asked,
 *            or the image of voxels if the transpose projector is asked.
 *
 *      in - input array (must be declared first, with right dimensions):
 *           it can be either the image of voxels if the direct projector is asked,
 *           or the data sinogram if the transpose projector is asked.
 *
 *      nx - size of the voxels image in x-direction.
 *
 *      ny - size of the voxels image in y-direction.
 *
 *      nz - size of the voxels image in z-direction.
 *
 *      s_scl - sampling step of voxels in each direction.
 *
 *      s_deg - B-spline degree of the image representation (ONLY DEGREE 3 IMPLEMENTED).
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      nv - detector size in v-direction (same as y-direction at initial position).
 *
 *      u_scl - sampling step of detector pixels in u-direction.
 *
 *      v_scl - sampling step of detector pixels in v-direction.
 *
 *      xos - coefficients of the 3D transform matrix, which "rotates" the
 *            voxels coordinates in the detector frame.
 *
 *      xsd - coefficients of the 3D transform matrix corresponding to the
 *            source position relative to the detector plane.
 *
 *      job - defines if the direct or transpose projector has to be applied
 *            (0 for direct / 1 for tranpose).
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <trtt_CB3D>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
extern int trtt_PB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job);

/*
 * Function: trtt_DD_PB3D
 *
 *      3D projector based on the Distance Driven model, in parallel beam geometry.
 *
 * Description:
 *
 *      The description of this function is the same as for <trtt_PB3D>.
 *
 * Parameters:
 *
 *      out - output array (must be declared first, with right dimensions):
 *            it can be either the data sinogram if the direct projector is asked,
 *            or the image of voxels if the transpose projector is asked.
 *
 *      in - input array (must be declared first, with right dimensions):
 *           it can be either the image of voxels if the direct projector is asked,
 *           or the data sinogram if the transpose projector is asked.
 *
 *      nx - size of the voxels image in x-direction.
 *
 *      ny - size of the voxels image in y-direction.
 *
 *      nz - size of the voxels image in z-direction.
 *
 *      s_scl - sampling step of voxels in each direction.
 *
 *      s_deg - B-spline degree of the image representation (ONLY DEGREE 3 IMPLEMENTED).
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      nv - detector size in v-direction (same as y-direction at initial position).
 *
 *      u_scl - sampling step of detector pixels in u-direction.
 *
 *      v_scl - sampling step of detector pixels in v-direction.
 *
 *      xos - coefficients of the 3D transform matrix, which "rotates" the
 *            voxels coordinates in the detector frame.
 *
 *      xsd - coefficients of the 3D transform matrix corresponding to the
 *            source position relative to the detector plane.
 *
 *      job - defines if the direct or transpose projector has to be applied
 *            (0 for direct / 1 for tranpose).
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <trtt_PB3D>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
extern int trtt_DD_PB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job);

/*
 * Function: trtt_CB3D
 *
 *      3D projector based on the B-spline model, in cone beam geometry.
 *
 * Description:
 *
 *      The description of this function is the same as for <trtt_PB3D>.
 *
 * Parameters:
 *
 *      out - output array (must be declared first, with right dimensions): it
 *            can be either the data sinogram if the direct projector is
 *            asked, or the image of voxels if the transpose projector is
 *            asked (keyword <job>).
 *
 *      in - input array (must be declared first, with right dimensions): it
 *           can be either the image of voxels if the direct projector is
 *           asked, or the data sinogram if the transpose projector is asked
 *           (keyword <job>).
 *
 *      nx - size of the voxels image in x-direction.
 *
 *      ny - size of the voxels image in y-direction.
 *
 *      nz - size of the voxels image in z-direction.
 *
 *      s_scl - sampling step of voxels in each direction.
 *
 *      s_deg - B-spline degree of the image representation (ONLY DEGREE 3
 *              IMPLEMENTED).
 *
 *      nu - detector size in u-direction (equivalent to z-direction).
 *
 *      nv - detector size in v-direction (same as y-direction at initial
 *           position).
 *
 *      u_scl - sampling step of detector pixels in u-direction.
 *
 *      v_scl - sampling step of detector pixels in v-direction.
 *
 *      xos - coefficients of the 3D transform matrix, which "rotates" the
 *            voxels coordinates in the detector frame.
 *
 *      xsd - coefficients of the 3D transform matrix corresponding to the
 *            source position relative to the detector plane.
 *
 *      job - defines if the direct or transpose projector has to be applied
 *            (0 for direct / 1 for tranpose).
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *
 * See also:
 *
 *      <trtt_PB3D>.
 *
 * MiTiV pages:
 *      <http://mitiv.univ-lyon1.fr/index.php?n=Tomographie.ModeleBspline3D>
 *
 */
extern int trtt_CB3D(double *out, double *in, long nx, long ny, long nz, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nu, long nv, double u_scl, double v_scl, double *xos, double *xsd, int job);

TRTT_END_C_DECLS

#endif /* TRTT_3D_PROJECTORS_H */

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
