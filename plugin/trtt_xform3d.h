/*
 * trtt_xform3d.h --
 *
 * Definitions for 3D transformations.
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

#ifndef TRTT_XFORM3D_H
#define TRTT_XFORM3D_H 1

#include <stdio.h>
#include "trtt_common.h"

TRTT_BEGIN_C_DECLS

/*
 * Title: TR-TT 3D Coordinate Transforms
 *
 *   Version 1.0
 *
 *   Copyright (C) 2011, the MiTiV project.
 *
 *   This package gives functions to manipulate 3D coordinate transforms.
 *   The definitions are available by:
 *   : #include <xform3d.h>
 */

/*
 * Typedef: xform3d_t
 *   Opaque structure used to store a 3D coordinate transform.
 *
 */
typedef struct _trtt_xform3d trtt_xform3d_t;

/*
 * Function: trtt_xform3d_new
 *   Allocate new instance of a 3D coordinate transform.
 *
 * Description:
 *   This function allocates a new instance of 3D coordinate transform
 *   and initializes it with the identity.
 *
 * Return:
 *   The address of a new transform instance or NULL in case
 *   of error.
 *
 * See:
 *   <trtt_xform3d_delete()>.
 */
extern trtt_xform3d_t *trtt_xform3d_new();

extern int trtt_xform3d_init(trtt_xform3d_t *xform);

/*
 * Function: trtt_xform3d_delete
 *   Destroys an instance of a 3D coordinate transform.
 *
 * Description:
 *   This function releases whatever ressources are allocated
 *   to store a 3D coordinate transform.
 *
 * Parameters:
 *   xform - the instance to delete.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern void trtt_xform3d_delete(trtt_xform3d_t *xform);


/*
 * Function: trtt_xform3d_reset
 *   Reset a 3D coordinate transform to the identity transform.
 *
 * Description:
 *   This function stores the coefficients of the identity in
 *   _xform_, an instance of a 3D coordinate transform.
 *
 * Parameters:
 *   xform - the coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern int trtt_xform3d_reset(trtt_xform3d_t *xform);

/*
 * Function: trtt_xform3d_copy
 *   Copy a 3D coordinate transform.
 *
 * Description:
 *   This function copies the parameters of a 3D coordinate transform
 *   from a source instance _src_ to a destination instance _dst_.
 *
 * Parameters:
 *   dst   - the destination coordinate transform;
 *   src   - the source coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern int trtt_xform3d_copy(trtt_xform3d_t *dst,
                             const trtt_xform3d_t *src);

/*
 * Function: trtt_xform3d_apply
 *   Apply a 3D coordinate transform to a vector of coordinates.
 *
 * Description:
 *   This function transforms a triplet of coordinates provided by a 3 element
 *   source vector {x,y,z} and stores the result in a destination vector.  The
 *   operation can be done "in-place"; that is, with _dst_ = _src_.
 *
 * Parameters:
 *   dst   - the destination coordinates as a 3-element vector {x',y',z'};
 *   xform - the coordinate transform;
 *   src   - the source coordinates as a 3-element vector {x,y,z}.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern int trtt_xform3d_apply(double dst[3],
                              const trtt_xform3d_t *xform,
                              const double src[3]);

/*
 * Function: trtt_xform3d_combine
 *   Combine two 3D coordinate transforms.
 *
 * Description:
 *   This function combines the effect of two 3D coordinate transforms, _left_
 *   and _right_, to form a transform whose effects are the same as applying
 *   the right transform followed by the left one.  The resulting transform is
 *   stored in the destination _dst_.  The operation can be done "in-place";
 *   that is, with _dst_ = _left_ and/or _dst_ = _right_.
 *
 * Parameters:
 *   dst   - the destination coordinate transform;
 *   left  - the coordinate transform of the left operand;
 *   right - the coordinate transform of the right operand.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern int trtt_xform3d_combine(trtt_xform3d_t *dst,
                                const trtt_xform3d_t *left,
                                const trtt_xform3d_t *right);

/*
 * Function: trtt_xform3d_invert
 *   Invert a 3D coordinate transform.
 *
 * Description:
 *   This function inverts the effect of the source 3D coordinate transform,
 *   _src_, and stores the result in the destination _dst_.  The operation can
 *   be done "in-place"; that is, with _dst_ = _src_.
 *
 * Parameters:
 *   dst   - the destination coordinate transform;
 *   src   - the source coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern int trtt_xform3d_invert(trtt_xform3d_t *dst,
                               const trtt_xform3d_t *src);

/*
 * Function: trtt_xform3d_move
 *   Add a translation to a 3D coordinate transform.
 *
 * Description:
 *   This function modifies a coordinate transform A (given by _src_) to add a
 *   translation T (specified by the arguments _tx_, _ty_ and _tz_) to produce
 *   a resulting transform B (stored in _dst_).  The translation is applied
 *   *after* whatever coordinate transform is registered in _src_.  Formally
 *   this corresponds to
 *   :  B <-- T.A
 *   where T is the translation, A is the initial transform and B is the
 *   resulting transform.  The operation can be done "in-place"; that is, with
 *   _dst_ = _src_.
 * 
 * Parameters:
 *   dst - the destination transform;
 *   tx  - the translation along X-axis;
 *   ty  - the translation along Y-axis;
 *   tz  - the translation along Z-axis;
 *   src - the initial coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>.
 */
extern int trtt_xform3d_move(trtt_xform3d_t *dst,
                             double tx,
                             double ty,
                             double tz,
                             const trtt_xform3d_t *src);

/*
 * Function: trtt_xform3d_rotate
 *   Add a rotation to a 3D coordinate transform.
 *
 * Description:
 *   This function modifies a coordinate transform A (given by _src_) to add a
 *   rotation R (specified by the "quaternion" _qw_, _qx_, _qy_ and _qz_) to
 *   produce a resulting transform B (stored in _dst_).  The rotation is
 *   applied *after* whatever coordinate transform is registered in _src_.
 *   Formally this corresponds to
 *   :  B <-- R.A
 *   where R is the rotation, A is the initial transform and B is the
 *   resulting transform.  The operation can be done "in-place"; that is, with
 *   _dst_ = _src_.  The quaternion does not need to be normalized (use
 *   <trtt_xform3d_scale()> if you want to scale the coordinates).  If the
 *   center of the rotation is not (0,0,0), you must "move" the coordinate
 *   transform backward and then forward with <trtt_xform3d_move()>.
 *
 *    More informations on quaternions on
 *    <http://en.wikipedia.org/wiki/Quaternion>
 *
 * Parameters:
 *   dst - the destination transform;
 *   qw  - the 1st quaternion coefficient;
 *   qx  - the 2nd quaternion coefficient;
 *   qy  - the 3rd quaternion coefficient;
 *   qz  - the 4th quaternion coefficient;
 *   src - the initial coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>, <trtt_xform3d_move()>, <trtt_xform3d_scale()>.
 */
extern int trtt_xform3d_rotate(trtt_xform3d_t *dst,
                               double qw,
                               double qx,
                               double qy,
                               double qz,
                               const trtt_xform3d_t *src);

/*
 * Function: trtt_xform3d_scale
 *   Change the scale in a 3D coordinate transform.
 *
 * Description:
 *   This function modifies a coordinate transform A (given by _src_) to add a
 *   scaling transform S (specified by the argument _scl_) to produce
 *   a resulting transform B (stored in _dst_).  The scaling is applied
 *   *after* whatever coordinate transform is registered in _src_.  Formally
 *   this corresponds to
 *   :  B <-- S.A
 *   where S is the scaling transform, A is the initial transform and B is the
 *   resulting transform.  The operation can be done "in-place"; that is, with
 *   _dst_ = _src_.
 * 
 * Parameters:
 *   dst                 - the destination transform;
 *   x_scl, y_scl, z_scl - the scale factors on each direction;
 *   src                 - the initial coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>, <trtt_xform3d_move()>, <trtt_xform3d_rotate()>.
 */
extern int trtt_xform3d_scale(trtt_xform3d_t *dst,
                              double x_scl,
                              double y_scl,
                              double z_scl,
                              const trtt_xform3d_t *src);

/*
 * Function: trtt_xform3d_get_coefficients
 *   Get the coefficients of a 3D coordinate transform.
 *
 * Description:
 *   The 12 coefficients, a[0] ... a[11], correspond to the following 3D
 *   coordinate transform:
 *   :  x' = a[0]*x + a[1]*y +  a[2]*z +  a[3];
 *   :  y' = a[4]*x + a[5]*y +  a[6]*z +  a[7];
 *   :  z' = a[8]*x + a[9]*y + a[10]*z + a[11];
 *   where {x,y,z} and {x',y',z'} are respectively the coordinates before and
 *   after the transform.
 *
 * Parameters:
 *   dst - the destination array of 12 coefficients;
 *   src - the source coordinate transform.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>, <trtt_xform3d_set_coefficients()>.
 */
extern int trtt_xform3d_get_coefficients(double dst[12],
                                         const trtt_xform3d_t *src);

/*
 * Function: trtt_xform3d_set_coefficients
 *   Set the coefficients of a 3D coordinate transform.
 *
 * Description:
 *   See <trtt_xform3d_get_coefficients()> for the meaning of the
 *   coefficients.
 *
 * Parameters:
 *   dst - the destination coordinate transform;
 *   src - the source array of 12 coefficients.
 *
 * Returns:
 *
 *      A standard result: <TRTT_SUCCESS> on success; <TRTT_FAILURE> on error.
 *      In case of error, the global variable _errno_ is set with a relevant
 *      code.
 *
 * See:
 *   <trtt_xform3d_new()>, <trtt_xform3d_set_coefficients()>.
 */
extern int trtt_xform3d_set_coefficients(trtt_xform3d_t *dst,
                                         const double src[12]);

/*
 * Function: trtt_xform3d_print
 *   Print the matrix of a 3D coordinate transform.
 *
 * Description:
 *   This function prints the coefficients of a 3D coordinate transform
 *   in a human readable form.
 *
 * Parameters:
 *   output - the output file;
 *   xform  - the coordinate transform.
 *
 * See:
 *   <trtt_xform3d_new()>, <trtt_xform3d_get_coefficients()>.
 */
extern void trtt_xform3d_print(FILE *output, const trtt_xform3d_t *xform);

/*---------------------------------------------------------------------------*/

TRTT_END_C_DECLS

#endif /* TRTT_XFORM3D_H */

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
