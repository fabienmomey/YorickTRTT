/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_tomobj.i --
 *
 * 2D/3D Tomographic object package.
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
 * with a limited warran
ty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean that
 * it is complicated to manipulate, and that also therefore means that it is
 * reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more generally,
 * to use and operate it in the same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

/* REQUIREMENTS ============================================================== */

extern trtt_tomobj_info;
/* DOCUMENT trtt_tomobj_info

   Ytrtt_tomobj.i is the module of TRTT code for creating and manipulating
   tomographic objects. The concept is to have an "object" which carries all the
   parameters (dimensions, sampling steps) of the object of interest, including
   the image of samples values. As we use a b-spline representation, the object
   also stores spline interpolation coefficients to pass from the spline
   coefficients to the samples of the image. It is stored as a hashtab with the
   following entries:

   OBJECT
   |
   |_____ voxels                                     => image of samples values.   
   |
   |_____ nx, ny, nz, nt                             => dimensions.
   |
   |_____ s_scl, t_scl                               => sampling steps.
   |
   |_____ xi (XFORM3D)                               => scaling transform.
   |
   |_____ wx, wy, wz, wt, (spl_interp_coefs)         => interpolation coefficients.
   |
   |_____ s_deg, t_deg                               => spline degrees.
   |
   |_____ footprint, cum_footprint, size_footprint,
          step_footprint, coord_footprint            => footprint parameters
                                                        (monodimensionnal profile
                                                        of the projection
                                                        of a b-spline).

   The object of interest is generalized as a 4D object : spatio-temporal image
   (a sequence of NT 3D images NX x NY x NZ), so that we can make dynamic
   tomography. It is the most general case, and then it is also allowed to
   create static object, of several dimensions (2D or 3D). XI makes pass from
   the voxel coordinate system of the object to the metric coordinate
   system. Starting from this type of tomographic OBJECT, associated with one
   DATA object, we are able to build a projector. The standard model will be our
   Separable b-spline projector (Distance Driven will be also availbale), which
   needs a FOOTPRINT to be integrated on the detector plane. This FOOTPRINT is
   stored in the TOMOGRAPHIC OBJECT, as a monodimensionnal profile (due to the
   separability of the footprint on the detector). It is a normalized b-spline
   (Unit integrale, unit voxel sampling step) of spatial degree S_DEG (T_DEG
   corresponds to the temporal spline degree). The projector will be in charge
   to scale the footprint after having projected its position on the
   detector. Other entries are associated with the FOOTPRINT:
   
   - its number of values SIZE_FOOTPRINT ;
   
   - its centered normalized coordinates COORD_FOOTPRINT ;
   
   - the sampling step of these normalized coordinates STEP_FOOTPRINT ;
   
   - its cumulative sum CUM_FOOTPRINT (for integration on pixels detector).
   
   A TOMOGRAPHIC OBJECT is created with the function "trtt_tomobj_new". The
   arguments passed to the function are used to initialize all entries of such
   object structure.  XI is initialized to position the origin of the metric
   coordinate system at the center of the object (we already have all parameters
   to do that). The interpolation coefficients WX, WY, WZ and WT are calculated
   according to the specified spline degrees. All FOOTPRINT entries are
   calculated.

   Then the TOMOGRAPHIC OBJECT can be copied, saved and loaded from a YHD file
   (see doc of yeti_yhdf.i). Then new parameters can be set to the current
   object to change its features. We get the entries the same way as to get the
   keys of a hashtab.

   List of functions of the module (see the help of these functions):

   - trtt_tomobj_new()
   - trtt_tomobj_copy()
   - trtt_tomobj_save()
   - trtt_tomobj_load()
   - trtt_tomobj_set_voxels()
   - trtt_tomobj_set_temporal_slice()
   - trtt_tomobj_set_scl()
   - trtt_tomobj_set_spline_deg()
   - trtt_tomobj_set_size_footprint()
   - trtt_tomobj_get_vox_coordinates()
   - trtt_tomobj_is_obj()
   

   SEE ALSO: trtt_tomobj_new, trtt_tomdata_info.
*/

/* GLOBALS =================================================================== */

AXIS_X = 0n;
AXIS_Y = 1n;
AXIS_Z = 2n;
AXIS_T = 3n;

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

func _monodim_bspline(m, deg, status)
/* DOCUMENT _monodim_bspline, m, deg, status

   Return a monodimensionnal b-spline profile of degree DEG. The result is given
   as a hashtab containing the array of M values of the b-spline, the associated
   coordinates and the sampling step. The b-spline profile is normalized, that
   is to say that its integrale and the voxel sampling ara unity.  The b-spline
   profile is created with the facilities of the SPL library and the b-spline
   interpolation tools. STATUS is the error controler.

   SEE ALSO spl_interp_coefs_new_spline.
 */
{
    if(trtt_error_query(status)) {
        return TRTT_FAILURE;
    }

    n = 2*deg+3;

    /* Create interpolation coordinates and normalize them */
    x = indgen(m)-1;
    step = double(n-1)/double(m-1);
    /* Scaling normalized coordinates */
    x *= step;

    /* Create the nodes */
    c = array(double, n);
    midpoint = long(floor(double(n+1)/2.0));
    c(midpoint) = 1.0;
    /* The middle node is set to 1 and the others to 0 => there is a single
       Bspline for interpolation on a fine grid => the interpolation will "draw"
       this Bspline */

    /* Create the Bspline interpolator */
    ws = spl_interp_coefs_new_spline(x,n,deg);

    /* Centering coordinates */
    x -= 0.5*double(n-1);
    bspline = ws(c);

    return h_new(bspline=bspline, coord=x, step=step);
}

func _cumsum(y, m, step)
/* DOCUMENT _cumsum, y, m, step

   Return the cumulative sum of the array Y of M samples with sampling step
   STEP.
 */
{
    ycum = array(double, m);
    ycum(1) = 0.0;
    for (i=1; i<=m; ++i) {
        ycum(i) = ycum(i-1)+0.5*step*(y(i-1)+y(i));
    }
    return ycum;
}

func _tomobj_make_xi(obj, status)
/* DOCUMENT _tomobj_make_xi, obj, status

   Create and initialize XI transform of the tomographic object OBJ. XI makes
   pass from the voxel coordinate system to the metric coordinate system, the
   origin of which is taken at the center of the object image. STATUS is the
   error controler.

   SEE ALSO: trtt_xform3d_new.
 */
{
    if(trtt_error_query(status)) {
        return TRTT_FAILURE;
    }
    
    /* Get dimensions and scaling parameters */
    nx = obj.nx;
    ny = obj.ny;
    nz = obj.nz;
    s_scl = obj.s_scl;

    /* Calculate the centering parameters */
    /* FIXME: CENTERS OF BLOBS ARE TAKEN INTO ACCOUNT */
    tx = -0.5*(nx-1); //-0.5*nx;
    ty = -0.5*(ny-1); //-0.5*ny;
    tz = -0.5*(nz-1); //-0.5*nz;

    /* Allocate xi transform */
    xi = trtt_xform3d_new();

    /* Add the centering translation */
    trtt_xform3d_move, xi, tx, ty, tz, xi;
    /* Scale the xi transform with sampling parameters */
    trtt_xform3d_scale, xi, s_scl, s_scl, s_scl, xi;
    
    h_set, obj, xi=xi;

    return TRTT_SUCCESS;
}

func _make_interpolator(obj, axis, status)
/* DOCUMENT _make_interpolator, obj, axis, status

   Create the b-spline interpolation coefficients to pass from the b-spline
   coefficients C_k to the voxel samples F_k of the tomographic object OBJ. AXIS
   indicates for which axis direction the interpolator is created : x, y, z or t
   axis. The coordinates on which the values are interpolated are the same as
   the coordinates of the b-spline coefficients, so it simply consists in the
   convolution of the C_k coefficients by a b-spline profile. The degree DEG of
   the b-spline profile and the number interpolation coordinates are extracted
   from OBJ (S_DEG, T_DEG, NX, NY, NZ, NT). The interpolator is created with the
   facilities of the SPL library and the b-spline interpolation tools. STATUS
   the error controler.

   SEE ALSO: trtt_tomobj_info, spl_interp_coefs_new_spline.
 */
{
    if(trtt_error_query(status)) {
        return TRTT_FAILURE;
    }

    if (axis == AXIS_X) {
        n = obj.nx;
        deg = obj.s_deg;
    }
    else if (axis == AXIS_Y) {
        n = obj.ny;
        deg = obj.s_deg;
    }
    else if (axis == AXIS_Z) { 
        n = obj.nz;
        deg = obj.s_deg;
    }
    else if (axis == AXIS_T) { 
        n = obj.nt;
        deg = obj.t_deg;
    }

    /* Create new interpolator */
    w = spl_interp_coefs_new_spline(indgen(n)-1.0, n, deg);

    /* Store new interpolator */
    if (axis == AXIS_X) {
        h_set, obj, wx=w;
    }
    else if (axis == AXIS_Y) {
        h_set, obj, wy=w;
    }
    else if (axis == AXIS_Z) { 
        h_set, obj, wz=w;
    }
    else if (axis == AXIS_T) { 
        h_set, obj, wt=w;
    }

    return TRTT_SUCCESS;
}

func _make_footprint(obj, status)
/* DOCUMENT _make_footprint, obj, status

   Create the FOOTPRINT of the tomographic object OBJ, and all the associated
   entries: CUM_FOOTPRINT, COORD_FOOTPRINT, STEP_FOOTPRINT. The FOOTPRINT is a
   monodimensionnal profile (due to the separability of the footprint on the
   detector). It is a normalized b-spline (Unit integrale, unit voxel sampling
   step) of spatial degree S_DEG extracted from OBJ. STATUS the error controler.

   SEE ALSO: trtt_tomobj_info, _monodim_bspline.
 */
{
    if(trtt_error_query(status)) {
        return TRTT_FAILURE;
    }

    /* Get footprint parameters*/
    size_footprint = obj.size_footprint;
    s_deg = obj.s_deg;
    

    /* Calculate values  */
    footprint_features = _monodim_bspline(size_footprint, s_deg, status);
    if (!is_hash(footprint_features) && footprint_features==TRTT_FAILURE) {
        return TRTT_FAILURE;
    }
    footprint = footprint_features.bspline;
    coord_footprint = footprint_features.coord;
    step_footprint = footprint_features.step;
    

    /* Calculate partial sums */
    cum_footprint = _cumsum(footprint, size_footprint, step_footprint);

    /* Set footprint parameters */
    h_set, obj, footprint = footprint;
    h_set, obj, cum_footprint = cum_footprint;
    h_set, obj, step_footprint = step_footprint;
    h_set, obj, coord_footprint = coord_footprint;

    return TRTT_SUCCESS;
}

/* USER FUNCTIONS ============================================================ */

func trtt_tomobj_new(voxels, s_scl, t_scl, s_deg, t_deg, size_footprint, status=)
/* DOCUMENT trtt_tomobj_new, voxels, s_scl, t_scl, s_deg, t_deg, size_footprint, status=

   Create and initialize a TOMOGRAPHIC OBJECT, according to the following parameters:
   
   - VOXELS:         image of the object of interest.
   - S_SCL:          spatial sampling step.
   - T_SCL:          temporal sampling step. 
   - S_DEG:          spatial spline degree.
   - T_DEG:          temporal spline degree.
   - SIZE_FOOTPRINT: number of stored samples of the footprint.

   The image VOXELS can be a 3-dimensional array (for 2D tomography => 2 +1) or
   a 4-dimensional array (for 3D tomography => 3 +1). The last additional
   dimension corresponds to the temporal dimension and is mandatory even in the
   static tomography case for which it must be set to 1.

   A STATUS error controler can be given (see trtt_error_info). If not it is
   created in the function and if an error occurs, it will be directly
   displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_copy.
 */
{
    extern AXIS_X;
    extern AXIS_Y;
    extern AXIS_Z;
    extern AXIS_T;
    
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) { /* Check if there is already an error */
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init(void);
    }

    /* Check pararmeters */
    // Check voxels
    if (!is_voxels(voxels)) {
        trtt_error_set, status, TRTT_ERR_NOT_VOXELS, disp=flag_disp_status;
        return;
    }
    // Check scales
    if (!is_tomobj_scales(s_scl,t_scl)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return;
    }
    // Check spline degrees
    if (!is_tomobj_degs(s_deg,t_deg)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return;
    }
    // Check footprint size
    if (!is_1dim(size_footprint)) {
        trtt_error_set, status, TRTT_ERR_BAD_DIM_VALUE, disp=flag_disp_status;
        return;
    }

    /* Create tomographic object hash tab and fill the first parameters */
    dims = dimsof(voxels);
    // FIXME: We allow voxels of dim 2 (2D static), dim 3 (3D static if t_scl=0), dim 3 (2D+t) or 4 (3D+t)
    // But nz and nt are still specified in the OBJECT
    if (dims(1) == 2) {
        nz = 1;
        nt = 1;
        nx = dims(2);
        ny = dims(3); 
    } else if (dims(1) == 3 & t_scl!=0.0) {
        nz = 1;
        nx = dims(2);
        ny = dims(3);
        nt = dims(4);
    } else if (dims(1) == 3 & t_scl==0.0) {
        nt = 1;
        nx = dims(2);
        ny = dims(3);
        nz = dims(4);
    } else if (dims(1) == 4) {
        nx = dims(2);
        ny = dims(3);
        nz = dims(4);
        nt = dims(5);  
    } else {
        trtt_error_set, status, TRTT_ERR_NOT_VOXELS, disp=flag_disp_status;
        return;
    }
    obj = h_new(voxels = voxels,
                nx = nx, ny = ny, nz = nz, nt = nt,
                s_scl = s_scl, t_scl = t_scl,
                s_deg = s_deg, t_deg = t_deg,
                size_footprint = size_footprint
        );

    /* Add xi transform */
    if (_tomobj_make_xi(obj, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* Add interpolators */
    if (_make_interpolator(obj, AXIS_X, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_Y, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_Z, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_T, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* Add footprint */
    if (_make_footprint(obj, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* Object not corrupted */
    // h_set, obj, corrupt_flag = TRTT_FALSE;

    return obj;    
}

func trtt_tomobj_copy(obj, status=)
/* DOCUMENT trtt_tomobj_copy, obj, status=

   Copy the tomographic object OBJ. Return a TOMOGRAPHIC OBJECT with the same
   entries as OBJ. A STATUS error controler can be given. If not it is created
   in the function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_tomobj_info.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    return h_copy(obj, 1);
}

func trtt_tomobj_save(DSTFILE, obj, status=, overwrite=)
/* DOCUMENT trtt_tomobj_save, DSTFILE, obj, status=, overwrite=

   Save the tomographic object OBJ in a YHD file DSTFILE (see doc of yeti_yhdf.i
   to learn more about YHD files). A STATUS error controler can be given. If not
   it is created in the function and if an error occurs, it will be directly
   displayed. If keyword OVERWRITE is true and file DSTFILE already exists, the
   new file will (silently) overwrite the old one; othwerwise, file DSFILE must
   not already exist (default behaviour).

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_load, yhd_save.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* DSTFILE must be a string */
    if (!is_string(DSTFILE)) {
        trtt_error_set, status, TRTT_ERR_NOT_A_STRING, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }
    
    yhd_save, DSTFILE, obj, overwrite=overwrite;
    return TRTT_SUCCESS;
}

func trtt_tomobj_load(SRCFILE, status=)
/* DOCUMENT trtt_tomobj_load, SRCFILE, status=

   Load a TOMOGRAPHIC OBJECT from YHD file SRCFILE (see doc of yeti_yhdf.i to
   learn more about YHD files). A STATUS error controler can be given. If not it
   is created in the function and if an error occurs, it will be directly
   displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_save, yhd_restore.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* SRCFILE must be a string */
    if (!is_string(SRCFILE)) {
        trtt_error_set, status, TRTT_ERR_NOT_A_STRING, disp=flag_disp_status;
        return;
    }

    /* Restore hash tab from file */
    obj = yhd_restore(SRCFILE);

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    return obj;
}

func trtt_tomobj_set_voxels(obj, voxels, status=)
/* DOCUMENT trtt_tomobj_set_voxels, obj, voxels, status=

   Set a voxel image VOXELS in the tomographic object OBJ. If dimensions have
   changed, related entries in OBJ are modified. A STATUS error controler can be
   given. If not it is created in the function and if an error occurs, it will
   be directly displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_set_temporal_slice.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }
    
    /* Check voxels */
    if (!is_voxels(voxels)) {
        trtt_error_set, status, TRTT_ERR_NOT_VOXELS, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Set new voxels */
    h_set, obj, voxels=voxels;

    dims = dimsof(voxels);
    // FIXME: We allow voxels of dim 2 (2D static), dim 3 (3D static if t_scl=0), dim 3 (2D+t) or 4 (3D+t)
    // But nz and nt are still specified in the OBJECT
    if (dims(1) == 2) {
        nz = 1;
        nt = 1;
        nx = dims(2);
        ny = dims(3); 
    } else if (dims(1) == 3 & t_scl!=0.0) {
        nz = 1;
        nx = dims(2);
        ny = dims(3);
        nt = dims(4);
    } else if (dims(1) == 3 & t_scl==0.0) {
        nt = 1;
        nx = dims(2);
        ny = dims(3);
        nz = dims(4);
    } else if (dims(1) == 4) {
        nx = dims(2);
        ny = dims(3);
        nz = dims(4);
        nt = dims(5);  
    } else {
        trtt_error_set, status, TRTT_ERR_NOT_VOXELS, disp=flag_disp_status;
        return;
    }

    /* Set new dimensions */
    h_set, obj, nx=nx;
    h_set, obj, ny=ny;
    h_set, obj, nz=nz;
    h_set, obj, nt=nt;
    
    /* Add xi transform */
    if (_tomobj_make_xi(obj, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return TRTT_FAILURE;
    }

    /* Add interpolators */
    if (_make_interpolator(obj, AXIS_X, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_Y, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_Z, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_T, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return TRTT_FAILURE;
    }

    return TRTT_SUCCESS;
}

func trtt_tomobj_set_temporal_slice(obj, slice, t_index, status=)
/* DOCUMENT trtt_tomobj_set_temporal_slice, obj, slice, t_index, status=

   Set a temporal SLICE of a voxel image in the tomographic object
   OBJ. Dimensions have to be compatible with the spatial dimensions of OBJ. A
   STATUS error controler can be given. If not it is created in the function and
   if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_set_voxels.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }
    
    /* Check slice */
    nx = obj.nx;
    ny = obj.ny;
    nz = obj.nz;
    
    if (!is_slice(slice, nx, ny, nz)) {
        trtt_error_set, status, TRTT_ERR_NOT_VOXELS, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Check temporal index */
    if (!is_tomobj_t_index(t_index, obj.nt)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Copy slice */
    obj.voxels(..,t_index) = slice;

    return TRTT_SUCCESS;
}

func trtt_tomobj_set_scl(obj, s_scl, t_scl, status=)
/* DOCUMENT trtt_tomobj_set_scl, obj, s_scl, t_scl, status=

   Set new sampling steps S_SCL and T_SCL in the tomographic object OBJ. Related
   entries in OBJ such XI are modified. A STATUS error controler can be
   given. If not it is created in the function and if an error occurs, it will
   be directly displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_set_voxels.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Check dimensions */
    if (!is_tomobj_scales(s_scl,t_scl)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Set new scale parameters */
    h_set, obj, s_scl=s_scl;
    h_set, obj, t_scl=t_scl;

    /* Add xi transform */
    if (_tomobj_make_xi(obj, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return TRTT_FAILURE;
    }

    return TRTT_SUCCESS;
}

func trtt_tomobj_set_spline_deg(obj, s_deg, t_deg, status=)
/* DOCUMENT trtt_tomobj_set_spline_deg, obj, s_deg, t_deg, status=

   Set new spline degrees S_DEG and T_DEG in the tomographic object OBJ. Related
   entries in OBJ such FOOTPRINT and INTERPOLATORS are modified. A STATUS error
   controler can be given. If not it is created in the function and if an error
   occurs, it will be directly displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_set_scl.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Check dimensions */
    if (!is_tomobj_degs(s_deg,t_deg)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Set new scale parameters */
    h_set, obj, s_deg=s_deg;
    h_set, obj, t_deg=t_deg;

    /* Add interpolators */
    if (_make_interpolator(obj, AXIS_X, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_Y, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_Z, status) == TRTT_FAILURE ||
        _make_interpolator(obj, AXIS_T, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* Add footprint */
    if (_make_footprint(obj, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    return TRTT_SUCCESS;
}

func trtt_tomobj_set_size_footprint(obj, size_footprint, status=)
/* DOCUMENT trtt_tomobj_set_size_footprint, obj, size_footprint, status=

   Set new footprint size SIZE_FOOTPRINT in the tomographic object OBJ. Related
   entries in OBJ such FOOTPRINT are modified. A STATUS error controler can be
   given. If not it is created in the function and if an error occurs, it will
   be directly displayed.

   SEE ALSO: trtt_tomobj_info, trtt_tomobj_set_spline_deg.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* obj must be a tomobj */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Check dimensions */
    if (is_1dim(size_footprint) == TRTT_FAILURE) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Set new scale parameters */
    h_set, obj, size_footprint=size_footprint;
    
    /* Add footprint */
    if (_make_footprint(obj, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    return TRTT_SUCCESS;
}

func trtt_tomobj_get_vox_coordinates(obj, status=)
/* DOCUMENT trtt_tomobj_get_vox_coordinates, obj, status=

   Return a hashtab filled with the voxel coordinate in the metric space.

   COORD
   |
   |_ x
   |_ y
   |_ z

   A STATUS error controler can be given. If not it is created in the function
   and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_tomobj_info.
 */    
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        if(trtt_error_query(status)) {
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        status = trtt_error_init();
    }

    /* data must be a tomdata */
    if (!is_tomobj(obj)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    nx = obj.nx;
    ny = obj.ny;
    nz = obj.nz;
    s_scl = obj.s_scl;

    x = (indgen(nx)-0.5*(nx+1))*s_scl;
    y = (indgen(ny)-0.5*(ny+1))*s_scl;
    z = (indgen(nz)-0.5*(nz+1))*s_scl;

    return h_new(x=x,y=y,z=z);
}
trtt_vox_coord = trtt_tomobj_get_vox_coordinates;

func trtt_tomobj_is_obj(obj)
/* DOCUMENT trtt_tomobj_is_obj, obj

   Check if OBJ is a TOMOGRAPHIC OBJECT. Return true or false.

   SEE ALSO: trtt_tomdata_info.
 */
{
    list_param = ["voxels", "nx", "ny", "nz", "nt",
                  "s_scl", "t_scl",
                  "s_deg", "t_deg",
                  "size_footprint",
                  "xi",
                  "wx", "wy", "wz", "wt",
                  "footprint", "cum_footprint", "coord_footprint"];

    /* Is it a hash tab ?*/
    if (!is_hash(obj)) {
        return TRTT_FALSE;
    }

    /*FIXME: ADD A COMPARISON OF HASH TAB KEYWORDS*/
    
    return TRTT_TRUE;
}
is_tomobj = trtt_tomobj_is_obj; // Define an easier alias

/* INTERNAL LOW LEVEL FUNCTIONS ============================================= */

func _trtt_tomobj_is_dims(n1,n2,n3,n4)
/* DOCUMENT _trtt_tomobj_is_dims, n1, n2, n3, n4

   Check if N1, N2, N3 and N4 are dimensions.

   SEE ALSO: _trtt_is_1dim. 
*/
{
    if (!is_1dim(n1) ||
        !is_1dim(n2) ||
        !is_1dim(n3) ||
        is_1dim(n4) ) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_tomobj_dims = _trtt_tomobj_is_dims;


func _trtt_tomobj_is_scales(scl1, scl2)
/* DOCUMENT _trtt_tomobj_is_scales, scl1, scl2

   Check if SCL1 and SCL2 are sampling steps.

   SEE ALSO: _trtt_is_1scale. 
*/
{
    if (!is_1scale(scl1) ||
        !is_1scale(scl2)) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_tomobj_scales = _trtt_tomobj_is_scales;

func _trtt_tomobj_is_degs(deg1, deg2)
/* DOCUMENT _trtt_tomobj_is_degs, scl1, scl2

   Check if DEG1 and DEG2 are spline degrees.

   SEE ALSO: _trtt_is_1deg. 
*/
{
    if (!is_1deg(deg1) ||
        !is_1deg(deg2) ) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_tomobj_degs = _trtt_tomobj_is_degs;

func _trtt_tomobj_is_t_index(t_index, nt)
/* DOCUMENT _trtt_tomobj_is_t_index, t_index, nt

   Check if T_INDEX is a temporal index included into the right range [0 NT-1].
*/
{
    if (!is_scalar(t_index) || !is_integer(t_index) || !(t_index>=0 & t_index<=nt-1)) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_tomobj_t_index = _trtt_tomobj_is_t_index;

func _trtt_is_voxels(voxels)
/* DOCUMENT _trtt_is_voxels, voxels

   Check if VOXELS is an image of a tomographic object. It can be 3- or 4-
   dimensionnal, depending whether we do 2D or 3D tomography. The temporal
   dimension is mandatory even if we do static tomography (NT=1).
*/
{
    /* Check parameters */
    if (!is_array(voxels) || !is_real(voxels)) {
        return TRTT_FALSE;
    }
    dims = dimsof(voxels);
    // Check number of dimensions
    if (dims(1) != 2 & dims(1) != 3 & dims(1) != 4) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_voxels = _trtt_is_voxels;

func _trtt_is_slice(slice, nx, ny, nz)
/* DOCUMENT _trtt_is_voxels, voxels

   Check if SLICE is an image of a tomographic object. It can be 2- or 3-
   dimensionnal, depending whether we do 2D or 3D tomography. Then its
   dimensions have to be NX x NY or NX x NY x NZ.

*/
{
    /* Check parameters */
    if (!is_array(slice) || !is_real(slice)) {
        return TRTT_FALSE;
    }
    dims = dimsof(slice);
    // Check number of dimensions
    if (dims(1) != 2 & dims(1) != 3) {
        return TRTT_FALSE;
    }

    if ((dims(1) == 2) & (dims(2)!=nx || dims(3)!=ny)) return TRTT_FALSE;
    if ((dims(1) == 3) & (dims(2)!=nx || dims(3)!=ny || dims(4)!=nz)) return TRTT_FALSE;
        
    return TRTT_TRUE;
}
is_slice = _trtt_is_slice;


// _trtt_is_1deg = is_1deg;
// _trtt_is_1scale = is_1scale;
// _trtt_is_1dim = is_1dim;
