/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_tomdata.i --
 *
 * 1D/2D Tomographic data package.
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

extern trtt_tomdata_info;
/* DOCUMENT trtt_tomdata_info

   Ytrtt_tomdata.i is the module of TRTT code for creating and manipulating
   tomographic data. The concept is to have a DATA "object" which carries not
   only the data itself (that is to say the sinogram or projection), but also
   all the calibration system parameters : source position, detector dimensions
   and position, sampling steps. The calibration data is given as 3D spatial
   transforms (translation, rotation, scaling), which allow to pass from one
   coordinate system to another (for example from the source to the detector). A
   single DATA object corresponds to a single projection (1 source position). It
   is stored as a hashtab with the following entries:

   DATA
   |
   |_____ pixels                  => sinogram.      
   |
   |_____ nu, nv                  => dimensions of the detector.
   |
   |_____ u_scl, v_scl            => sampling steps. 
   |
   |_____ t_index                 => time coordinate index.
   |
   |_____ xi, xos, xotilt, xsd    => spatial transforms (calibration data).

   T_INDEX refers to the instant when the data is acquired from the tomographic
   object. It is used for dynamic tomography to link the projection to a
   particular slice of the sequence of images to reconstruct.
   
   The transforms XI, XOS, XOTILT, XSD, are XFORM3D objects, opaque structure
   storing coefficients of a spatial transform, developped in the SPL library
   (refer to it to know how to manipulate such objects).  
   - XI makes pass from the pixel coordinate system of the detector to the
     metric coordinate system of the detector plane. 
   - XOS make pass from the tomographic objet coordinate system to the source
     coordinate system.
   - XOTILT specify a tilt transform of the tomographic object at the instant
     T_INDEX of acquisition. It is defined in the DATA because it relies on a
     particular instant of the acquisition, and hence is related to the
     calibration of the system.
   - XSD make pass from the source coordinate system to the detector coordinate
     system. Its coefficients is used to project a voxel on the detector plane.

   A DATA object is created with the function "trtt_tomdata_new". The arguments
   passed to the function are used to initialize all entries of such object
   structure. The XFORM3D transforms are set to identity except for XI which is
   initialized to position the origin of the metric coordinate system at the
   center of the detector plane (we already have all parameters to do that).

   Then the DATA object can be copied, saved and loaded from a YHD file (see doc
   of yeti_yhdf.i). Then new parameters can be set to the current object to
   change its features. We get the entries the same way as to get the keys of a
   hashtab. For example we can get the transforms to add rotations, translations
   or scalings, for calibrating the projection.

   Here are some advice about how to modify the XFORM3D transforms to get a good
   calibration of the system, depending on the geometry of projection:

   - PARALLEL BEAM GEOMETRY:

   It is important to notice that the detector plane could not be orthogonal to
   the ray trajectory (=> distorsion factor for spline projector). Hence XSD
   could define a rotation of the detector coordinate system relative to the
   source coordinate system, to get the incidence of the rays on the
   detector. It is then combined with a translation in the detector axis
   directions to position the projected voxel points relative to the origin of
   the detector plane. The focale is technically useless (translation of the
   coordinate system in the rays direction), but you can spécify it to be sure.

   - CONE BEAM GEOMETRY:

   Due to the divergent trajectories of the rays, it is simplier that the source
   and detector coordinate systems are only translated one with respect to the
   other (=> it will simplifies the calculation of the distorsion
   factors). Hence XSD defines only a translation in the detector axis
   directions to position the projected voxel points relative to the origin of
   the detector plane, and a translation in the rays direction which defines the
   orthogonal projection of the source on the detector plane (focale).

   List of functions of the module (see the help of these functions):

   - trtt_tomdata_new()
   - trtt_tomdata_copy()
   - trtt_tomdata_save()
   - trtt_tomdata_load()
   - trtt_tomdata_set_pixels()
   - trtt_tomdata_set_scl()
   - trtt_tomdata_set_t_index()
   - trtt_tomdata_get_pix_coordinates()
   - trtt_tomdata_is_data()

   SEE ALSO: trtt_tomdata_new, trtt_tomobj_info.
*/

/* PRIVATE FUNCTIONS ==================================================== */

func _tomdata_make_xi(data, status)
/* DOCUMENT _tomdata_make_xi, data, status

   Create and initialize XI transform of DATA object. XI makes pass from the
   pixel detector coordinate system to the metric coordinate system, the origin
   of which is taken at the center of the detector plane. STATUS is the error
   controler.

   SEE ALSO: trtt_xform3d_new.
 */
{
    if(trtt_error_query(status)) {
        return TRTT_FAILURE;
    }
    
    /* Get dimensions and scaling parameters */
    nu = data.nu;
    nv = data.nv;
    u_scl = data.u_scl;
    v_scl = data.v_scl;
    // u_scl_inv = 1.0/u_scl;
    // v_scl_inv = 1.0/v_scl;

    /* Calculate the centering parameters */
    /* FIXME: CENTERS OF BLOBS ARE TAKEN INTO ACCOUNT */
    tu = -0.5*nu; //-0.5*(nu-1);
    tv = -0.5*nv; //-0.5*(nv-1);

    /* Allocate xi transform */
    xi = trtt_xform3d_new();

    /* Scale the xi transform with sampling parameters */
    // trtt_xform3d_scale, xi, u_scl_inv, v_scl_inv, 1.0, xi;
    /* Add the centering translation */
    // trtt_xform3d_move, xi, tu, tv, 0.0, xi;
   
    /* Add the centering translation */
    trtt_xform3d_move, xi, 0.0, tv, tu, xi;
    /* Scale the xi transform with sampling parameters */
    trtt_xform3d_scale, xi, 1.0, v_scl, u_scl, xi;

    
    h_set, data, xi=xi;

    return TRTT_SUCCESS;
}

func _make_xforms(data, status)
/* DOCUMENT _tomdata_xforms, data, status

   Create XOS, XOTILT and XSD transforms of DATA object. STATUS is the error
   controler.
 */
{
    if(trtt_error_query(status)) {
        return TRTT_FAILURE;
    }

    /* Allocate transforms */
    xos = trtt_xform3d_new();
    xotilt = trtt_xform3d_new();
    xsd = trtt_xform3d_new();

    h_set, data, xos=xos;
    h_set, data, xotilt=xotilt;
    h_set, data, xsd=xsd;
}

/* USER FUNCTIONS ============================================================ */

func trtt_tomdata_new(pixels, u_scl, v_scl, t_index, status=)
/* DOCUMENT trtt_tomdata_new, pixels, u_scl, v_scl, t_index, status=

   Create and initialize a DATA object, according to the following parameters:
   
   - PIXELS:  projection (or sinogram) data array.
   - U_SCL:   sampling step in the axial direction.
   - V_SCL:   sampling step in the transaxial direction. 
   - T_INDEX: time coordinate index (must be between 0 and NT-1, where NT is the
              number of temporal slices of the tomographic object).

   A STATUS error controler can be given (see trtt_error_info). If not it is
   created in the function and if an error occurs, it will be directly
   displayed.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_copy.
 */
{
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

    /* Check parameters */
    // Check pixels
    if (!is_pixels(pixels)) {
        trtt_error_set, status, TRTT_ERR_NOT_PIXELS, disp=flag_disp_status;
        return;
    }
    // Check scales
    if (!is_tomdata_scales(u_scl,v_scl)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return;
    }
    // Check t_index
    if (!is_tomdata_t_index(t_index)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return;
    }

    /* Create data hash tab and fill the first parameters */
    dims = dimsof(pixels);
    // FIXME: We allow pixels of dim 1 (2D) or 2 (3D)
    // But nu is still specified in the DATA
    if (dims(1) == 1) {
        nu = 1;
        nv = dims(2);
    } else {
        nu = dims(2);
        nv = dims(3);
    }
    data = h_new(pixels = pixels,
                 nu = nu, nv = nv,
                 u_scl = u_scl, v_scl = v_scl,
                 t_index = t_index
        );
    
    /* Add xi transform */
    if (_tomdata_make_xi(data, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* Add xforms */
    if (_make_xforms(data, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }
    
    return data;
}

func trtt_tomdata_copy(data, status=)
/* DOCUMENT trtt_tomdata_copy, data, status=

   Copy the DATA object. Return a DATA object with the same entries as DATA. A
   STATUS error controler can be given. If not it is created in the function and
   if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_tomdata_info.
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

    /* data must be a tomdata */
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return;
    }

    return h_copy(data, 1);
}

func trtt_tomdata_save(DSTFILE, data, status=, overwrite=)
/* DOCUMENT trtt_tomdata_save, DSTFILE, data, status=, overwrite=

   Save the DATA object in a YHD file DSTFILE (see doc of yeti_yhdf.i to learn
   more about YHD files). A STATUS error controler can be given. If not it is
   created in the function and if an error occurs, it will be directly
   displayed. If keyword OVERWRITE is true and file DSTFILE already exists, the
   new file will (silently) overwrite the old one; othwerwise, file DSFILE must
   not already exist (default behaviour).

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_load, yhd_save.
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

    /* data must be a tomdata */
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return TRTT_FAILURE;
    }
    
    yhd_save, DSTFILE, data, overwrite=overwrite;
    return TRTT_SUCCESS;
}

func trtt_tomdata_load(SRCFILE, status=)
/* DOCUMENT trtt_tomdata_load, SRCFILE, status=

   Load a DATA object from YHD file SRCFILE (see doc of yeti_yhdf.i to learn
   more about YHD files). A STATUS error controler can be given. If not it is
   created in the function and if an error occurs, it will be directly
   displayed.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_save, yhd_restore.
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
    data = yhd_restore(SRCFILE);

    /* data must be a tomdata */
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return;
    }

    return data;
}

func trtt_tomdata_set_pixels(data, pixels, status=)
/* DOCUMENT trtt_tomdata_set_pixels, data, pixels, status=

   Set new singoram data array PIXELS in the DATA object. If dimensions have
   changed, related entries in DATA are modified. A STATUS error controler can
   be given. If not it is created in the function and if an error occurs, it
   will be directly displayed.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_set_scl, trtt_tomdata_set_t_index.
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
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return TRTT_FAILURE;
    }
    
    /* Check pixels */
    if (!is_pixels(pixels)) {
        trtt_error_set, status, TRTT_ERR_NOT_PIXELS, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Set new voxels */
    h_set, data, pixels=pixels;
    
    dims = dimsof(pixels);
    // FIXME: We allow pixels of dim 1 (2D) or 2 (3D)
    // But nu is still specified in the DATA
    if (dims(1) == 1) {
        nu = 1;
        nv = dims(2);
    } else {
        nu = dims(2);
        nv = dims(3);
    }

    /* Set new dimensions */
    h_set, data, nu=nu;
    h_set, data, nv=nv;

    /* Add xi transform */
    if (_tomdata_make_xi(data, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return TRTT_FAILURE;
    }
    
    return TRTT_SUCCESS;
}

func trtt_tomdata_set_scl(data, u_scl, v_scl, status=)
/* DOCUMENT trtt_tomdata_set_scl, data, u_scl, v_scl, status=

   Set new sampling steps U_SCL and V_SCL in the DATA object. Related entries in
   DATA such XI are modified. A STATUS error controler can be given. If not it
   is created in the function and if an error occurs, it will be directly
   displayed.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_set_pixels, trtt_tomdata_set_t_index.
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
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Check dimensions */
    if (!is_tomdata_scales(u_scl,v_scl)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Set new scale parameters */
    h_set, data, u_scl=u_scl;
    h_set, data, v_scl=v_scl;

    /* Add xi transform */
    if (_tomdata_make_xi(data, status) == TRTT_FAILURE) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return TRTT_FAILURE;
    }

    return TRTT_SUCCESS;
}

func trtt_tomdata_set_t_index(data, t_index, status=)
/* DOCUMENT trtt_tomdata_set_t_index, data, t_index, status=

   Set new time coordinate index T_INDEX in the DATA object. A STATUS error
   controler can be given. If not it is created in the function and if an error
   occurs, it will be directly displayed.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_set_pixels, trtt_tomdata_set_scl.
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
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    /* Check t_index */
    if (!is_tomdata_t_index(t_index)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_VALUE, disp=flag_disp_status;
        return;
    }

    /* Set new t_index */
    h_set, data, t_index = t_index;

    return TRTT_SUCCESS;
}

func trtt_tomdata_get_pix_coordinates(data, status=)
/* DOCUMENT trtt_tomdata_get_pix_coordinates, data, status=

   Return a hashtab filled with the pixels detector coordinate in the metric
   space.

   COORD
   |
   |_ u
   |_ v

   A STATUS can be given. If not it is created in the function
   and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_tomdata_info.
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
    if (!is_tomdata(data)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    nu = data.nu;
    nv = data.nv;
    u_scl = data.u_scl;
    v_scl = data.v_scl;

    u = (indgen(nu)-0.5*(nu+1))*u_scl;
    v = (indgen(nv)-0.5*(nv+1))*v_scl;

    return h_new(u=u,v=v);
}
trtt_pix_coord = trtt_tomdata_get_pix_coordinates;

func trtt_tomdata_is_data(data)
/* DOCUMENT trtt_tomdata_is_data, data

   Check if DATA is a DATA object. Return true or false.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_is_data_set.
 */
{
    list_param = ["pixels", "nu", "nv",
                  "u_scl", "v_scl",
                  "t_index",
                  "xi", "xos", "xotilt", "xsd"];

    /* Is it a hash tab ?*/
    if (!is_hash(data)) {
        return TRTT_FALSE;
    }

    /*FIXME: ADD A COMPARISON OF HASH TAB KEYWORDS*/
    
    return TRTT_TRUE;
}
is_tomdata = trtt_tomdata_is_data;

func trtt_tomdata_is_data_set(datas)
/* DOCUMENT trtt_tomdata_is_data_set, datas

   Check if DATAS is a set of DATA objects (hashtab filled with several DATA
   objects + a list of DATA objects names). Return true or false.

   SEE ALSO: trtt_tomdata_info, trtt_tomdata_is_data.
 */
{
    /* Is it a hash tab ?*/
    if (!is_hash(datas)) {
        return TRTT_FALSE;
    }

    /*FIXME: ADD A COMPARISON OF HASH TAB KEYWORDS*/
    Yk_list = datas.Yk_list;
    if (is_void(datas.Yk_list) || !is_array(Yk_list) || !is_string(Yk_list)) {
        return TRTT_FALSE;
    }

    Ndata = numberof(Yk_list);

    for (k=1; k<=Ndata; ++k) {
        Yk = h_get(Y,Yk_list(k));

        /* Yk must be a tomdata */
        if (!is_tomdata(Yk)) {
            return TRTT_FALSE;
        }
    }
    
    return TRTT_TRUE;
}
is_tomdata_set = trtt_tomdata_is_data_set;


/* INTERNAL LOW LEVEL FUNCTIONS ============================================= */

func _trtt_tomdata_is_dims(n1,n2)
/* DOCUMENT _trtt_tomdata_is_dims, n1, n2

   Check if N1 and N2 are dimensions.

   SEE ALSO: _trtt_is_1dim. 
*/
{
        if (!is_1dim(n1) ||
            !is_1dim(n2)) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_tomdata_dims = _trtt_tomdata_is_dims;

func _trtt_tomdata_is_scales(scl1,scl2)
/* DOCUMENT _trtt_tomdata_is_scales, scl1, scl2

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
is_tomdata_scales = _trtt_tomdata_is_scales;

func _trtt_tomdata_is_t_index(t_index)
/* DOCUMENT _trtt_tomdata_is_t_index, t_index

   Check if T_INDEX is a temporal index.
*/
{
    if (!is_scalar(t_index) || !is_real(t_index)) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_tomdata_t_index = _trtt_tomdata_is_t_index;

func _trtt_is_pixels(pixels)
/* DOCUMENT _trtt_is_pixels, pixels

   Check if PIXELS is a sinogram data array. It can be 1- or 2- dimensionnal,
   depending whether we do 2D or 3D tomography.  
*/
{
    /* Check parameters */
    if (!is_array(pixels) || !is_real(pixels)) {
        return TRTT_FALSE;
    }
    dims = dimsof(pixels);
    // Check number of dimensions
    if (dims(1) != 1 & dims(1) != 2) {
        return TRTT_FALSE;
    }
    return TRTT_TRUE;
}
is_pixels = _trtt_is_pixels;
