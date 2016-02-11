/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_3D_projectors.i --
 *
 * TRTT 3D Radon projectors.
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

extern trtt_3D_projectors_info;

/* GLOBALS =================================================================== */

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

/* USER FUNCTIONS ============================================================ */

/*|---------------|*/
/*| 3D PROJECTORS |*/
/*|---------------|*/

func trtt_single_3D_projector(X, Yk, mode, status=)
/* DOCUMENT trtt_single_2D_projector, X, Yk, mode, status=

   Create a 3D projector with the Separable B-spline model, according to the
   geometry of projection MODE (TRTT_PARALLEL_BEAM or TRTT_FAN(CONE)_BEAM). X is
   the TOMOGRAPHIC OBJECT and Yk is the DATA object. The projector consists in a
   h_evaluator carrying references to X and Yk, associated with the adapted
   projection function. A STATUS can be given. If not it is created in the
   function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_3D_projectors_info, _trtt_single_3D_projector, h_evaluator.
 */    
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    /* Yk must be a tomdata */
    if (!is_tomdata(Yk)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return;
    }

    /* Extract needed parameters */
    /* - X */
    obj_xi = X.xi;
    nt = X.nt;
    /* - Yk */
    data_xi = Yk.xi;
    xotilt = Yk.xotilt;
    xos = Yk.xos;
    xsd = Yk.xsd;
    t_index = Yk.t_index;
    /* xform combinations */
    xos_plus = trtt_xform3d_combine(xotilt, obj_xi);
    trtt_xform3d_combine, xos_plus, xos, xos_plus;
    xos_coefs = xos_plus.a;
    data_xi_inv = trtt_xform3d_invert(data_xi);
    xsd_plus = trtt_xform3d_combine(data_xi_inv, xsd);
    xsd_coefs = xsd_plus.a;

    /* Check t_index */
    if (!(t_index>=0 & t_index<=nt)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
        return;
    }

    /* Create the evaluator */
    p = h_new(X=X,
              Yk=Yk,
              xos_coefs = xos_coefs,
              xsd_coefs = xsd_coefs);

    if (mode == TRTT_PARALLEL_BEAM) {
        h_set, p, projector = symlink_to_name("trtt_PB3D");
    }
    else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
        h_set, p, projector = symlink_to_name("trtt_CB3D");
    }
    else {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    h_evaluator, p, "_trtt_single_3D_projector";

    return p;
}
A_single_3D = trtt_single_3D_projector;

func _trtt_single_3D_projector(p, x, job)
/* DOCUMENT _trtt_single_3D_projector, p, x, job

   Evaluator (direct and tranpose) of a 3D projector using the Separable
   B-spline model. According to the geometry of projection, a call to the suited
   C-implemented projector is done.

   SEE ALSO: trtt_single_3D_projector, h_evaluator.
 */ 
{
    nx = p.X.nx;
    ny = p.X.ny;
    nz = p.X.nz;
    s_scl = p.X.s_scl;
    s_deg = p.X.s_deg;
    local footprint; eq_nocopy, footprint, p.X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, p.X.cum_footprint;
    size_footprint = p.X.size_footprint;
    step_footprint = p.X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, p.X.coord_footprint;
    nu = p.Yk.nu;        nv = p.Yk.nv; 
    u_scl = p.Yk.u_scl;  v_scl = p.Yk.v_scl;
    local xos_coefs; eq_nocopy, xos_coefs, p.xos_coefs;
    local xsd_coefs; eq_nocopy, xsd_coefs, p.xsd_coefs;
    __projector = p.projector;

    if (!job) {
        job=0n;
        out = array(double, nu, nv);
    } else if (job==1) {
        out = array(double, nx, ny, nz);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }
    
    if(__projector(out, x, nx, ny, nz, s_scl, s_deg,
                   footprint, cum_footprint,
                   size_footprint, step_footprint, coord_footprint,
                   nu, nv, u_scl, v_scl,
                   xos_coefs, xsd_coefs, job) == TRTT_FAILURE) {
        trtt_error_display_msg, code=TRTT_ERR_PROJECTOR_FAILURE;
        return;
    } else {
        return out;
    }             
}

func trtt_whole_3D_projector(X, Y, mode, status=)
/* DOCUMENT trtt_whole_3D_projector, X, Y, mode, status=

   Create a set of 3D projectors based on the Separable B-spline model, from a
   TOMOGRAPHIC OBJECT X and a set of DATA objects Y, according to the geometry
   of projection MODE. A STATUS can be given. If not it is created in the
   function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_3D_projectors_info, trtt_single_3D_projector.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    // FIXME: Check if Y is a hash table containing TOMDATA objects
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return;
    }

    Yk_list = Y.Yk_list;
    Ndata = numberof(Yk_list);

    /* Create then hash table which contains the whole projectors */
    Cops = h_new();

    for (k=1; k<=Ndata; ++k) {
        Yk = h_get(Y, Yk_list(k));

        Ck = A_single_3D(X, Yk, mode, status=status);
        if(is_void(Ck)) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }

        h_set, Cops, trtt_get_key_name(k, "C"), Ck;
    }

    Ck_list = h_keys(Cops);  
    h_set, Cops, Ck_list = Ck_list(sort(Ck_list));
    
    return Cops;
}
A_whole_3D = trtt_whole_3D_projector;

/*|-----------------------|*/
/*| 3D TIME INTERPOLATORS |*/
/*|-----------------------|*/

func trtt_single_3D_time_interpolator(X, Yk, status=)
/* DOCUMENT
 */    
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    /* Yk must be a tomdata */
    if (!is_tomdata(Yk)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return;
    }

    /* Extract needed parameters */
    /* - X */
    nt = X.nt;
    nx = X.nx;
    ny = X.ny;
    nz = X.nz;
    t_scl = X.t_scl;
    t_deg = X.t_deg;
    /* - Yk */
    t_index = Yk.t_index;

    /* Check t_index */
    if (!(t_index>=0 & t_index<=nt)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOOBJ.", disp=flag_disp_status;
        return;
    }

    /* Create the interpolator */
    wt = spl_interp_coefs_new_spline(t_index, nt, t_deg);

    p = h_new(wt = wt,
              nx = nx,
              ny = ny,
              nz = nz,
              interpolator = symlink_to_name("_trtt_3D_time_interpolator"));

    h_evaluator, p, "_trtt_single_3D_time_interpolator";

    return p;
}
B_single_3D = trtt_single_3D_time_interpolator;

func _trtt_single_3D_time_interpolator(p, x, job)
/* DOCUMENT
 */    
{
    wt = p.wt;
    __interpolator = p.interpolator;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    m = wt.m;
    n = wt.n; // = nt
    p = wt.p;
    j = wt.j + 1;
    w = wt.w;

    
    if (!job) {
        job=0n;
        out = array(double, nx, ny, nz);
    } else if (job==1) {
        out = array(double, nx, ny, nz, n);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    __interpolator, out, x, m, n, p, j, w, job;
    
    return out;
}

func _trtt_3D_time_interpolator(&out, x, m, n, p, j, w, job)
{
    // Allocate array for the result.
    if (!job) {
        job=0n;
    } else if (job==1) {
    } else {
        error, "unsupported JOB";
    }

    if (!job) {
        out = (w(-,-,-,)*x(..,j))(..,sum);
    } else if (job==1) {
        out(..,j) = w(-,-,-,)*x(..,-:1:p);
    } else {
        error, "unsupported JOB";
    }

    return TRTT_SUCCESS;
}

/*|----------------------------|*/
/*| DISTANCE DRIVEN PROJECTORS |*/
/*|----------------------------|*/

// func trtt_single_distance_driven_3D_projector(X, Yk, mode, status=)
// /* DOCUMENT trtt_single_distance_driven_3D_projector, X, Yk, mode, status=

//    This function computes a linear operator of algebraic Radon transform,
//    performing the "distance driven" method of DeMan & Basu (2004). It is created
//    from a tomobj X and a tomdata Yk. It can be performed both in "parallel" and
//    "fan beam" MODE of projection. The "distance driven" method is a combination
//    of the "voxel" and "ray driven" methods, considering voxels and detectors,
//    not as points anymore, but with the thickness of their support. The interest
//    is to prevent the artifacts that are generated in the "pixel driven"
//    projection and in the "ray driven" backprojection. The projector consists in
//    a h_evaluator carrying references to X and Yk, associated with the adapted
//    projection function. A STATUS can be given. If not it is created in the
//    function and if an error occurs, it will be directly displayed.

//    SEE ALSO: trtt_3D_projectors_info, _trtt_single_distance_driven_3D_projector, h_evaluator.
//  */    
// {
//     if (!is_void(status)) {
//         flag_disp_status = TRTT_FALSE;
//         /* A status variable is given : in case of error, status is updated to
//          * be passed to upper callers */
//         if(trtt_error_query(status)) {
//             return;
//         }
//     } else {
//         flag_disp_status = TRTT_TRUE;
//         /* We initialize a status variable for internal functions */
//         status = trtt_error_init();
//     }

//     /* X must be a tomobj */
//     if (!is_tomobj(X)) {
//         trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
//         return;
//     }

//     /* Yk must be a tomdata */
//     if (!is_tomdata(Yk)) {
//         trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
//         return;
//     }

//     /* Extract needed parameters */
//     /* - X */
//     obj_xi = X.xi;
//     nt = X.nt;
//     /* - Yk */
//     data_xi = Yk.xi;
//     xotilt = Yk.xotilt;
//     xos = Yk.xos;
//     xsd = Yk.xsd;
//     t_index = Yk.t_index;
//     /* xform combinations */
//     xos_plus = trtt_xform3d_combine(xotilt, obj_xi);
//     trtt_xform3d_combine, xos_plus, xos, xos_plus;
//     xos_coefs = xos_plus.a;
//     data_xi_inv = trtt_xform3d_invert(data_xi);
//     xsd_plus = trtt_xform3d_combine(data_xi_inv, xsd);
//     xsd_coefs = xsd_plus.a;

//     /* Check t_index */
//     if (!(t_index>=0 & t_index<=(nt-1))) {
//         trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
//         return;
//     }

//     /* Create the evaluator */
//     p = h_new(X=X,
//               Yk=Yk,
//               xos_coefs = xos_coefs,
//               xsd_coefs = xsd_coefs);

//     if (mode == TRTT_PARALLEL_BEAM) {
//         h_set, p, projector = symlink_to_name("trtt_DD_PB3D")
//     }
//     else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
//         h_set, p, projector = symlink_to_name("trtt_DD_CB3D")
//     }
//     else {
//         trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
//         return;
//     }

//     h_evaluator, p, "_trtt_single_distance_driven_3D_projector";

//     return p;
// }
// D_single_3D = trtt_single_distance_driven_3D_projector;

// func _trtt_single_distance_driven_3D_projector(p, x, job)
// /* DOCUMENT _trtt_single_distance_driven_3D_projector, p, x, job

//    Evaluator of the projector (direct and transpose) based on the "distance
//    driven" method of DeMan & Basu (2002), in "parallel beam" mode of projection.

//    SEE ALSO: trtt_single_distance_driven_3D_projector, h_evaluator.
//  */ 
// {
//     nx = p.X.nx;
//     ny = p.X.ny;
//     nz = p.X.nz;
//     s_scl = p.X.s_scl;
//     // s_deg = p.X.s_deg;
//     // local footprint; eq_nocopy, footprint, p.X.footprint;
//     // local cum_footprint; eq_nocopy, cum_footprint, p.X.cum_footprint;
//     // size_footprint = p.X.size_footprint;
//     // step_footprint = p.X.step_footprint;
//     // local coord_footprint; eq_nocopy, coord_footprint, p.X.coord_footprint;
//     nu = p.Yk.nu;        nv = p.Yk.nv; 
//     u_scl = p.Yk.u_scl;  v_scl = p.Yk.v_scl;
//     local xos_coefs; eq_nocopy, xos_coefs, p.xos_coefs;
//     local xsd_coefs; eq_nocopy, xsd_coefs, p.xsd_coefs;
//     __projector = p.projector;

//     if (!job) {
//         job=0n;
//         out = array(double, nu, nv);
//     } else if (job==1) {
//         out = array(double, nx, ny, nz);
//     } else {
//         trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
//         return;
//     }

//     if(__projector(out, x, nx, ny, nz, s_scl,
//                    nu, nv, u_scl, v_scl,
//                    xos_coefs, xsd_coefs, job) == TRTT_FAILURE) {
//         trtt_error_display_msg, code=TRTT_ERR_PROJECTOR_FAILURE;
//         return;
//     } else {
//         return out;
//     }             
// }

// func trtt_whole_distance_driven_3D_projector(X, Y, mode, status=)
// /* DOCUMENT trtt_whole_distance_driven_3D_projector, X, Y, mode, status=

//    Create a set of 3D projectors based on the Distance Driven model, from a
//    TOMOGRAPHIC OBJECT X and a set of DATA objects Y, according to the geometry
//    of projection MODE. A STATUS can be given. If not it is created in the
//    function and if an error occurs, it will be directly displayed.

//    SEE ALSO: trtt_3D_projectors_info, trtt_single_distance_driven_3D_projector.
//  */
// {
//     if (!is_void(status)) {
//         flag_disp_status = TRTT_FALSE;
//         /* A status variable is given : in case of error, status is updated to
//          * be passed to upper callers */
//         if(trtt_error_query(status)) {
//             return;
//         }
//     } else {
//         flag_disp_status = TRTT_TRUE;
//         /* We initialize a status variable for internal functions */
//         status = trtt_error_init();
//     }

//     /* X must be a tomobj */
//     if (!is_tomobj(X)) {
//         trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
//         return;
//     }

//     // FIXME: Check if Y is a hash table containing TOMDATA objects
//     if (!is_tomdata_set(Y)) {
//         trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
//         return;
//     }

//     Yk_list = Y.Yk_list;
//     Ndata = numberof(Yk_list);

//     /* Create then hash table which contains the whole projectors */
//     Cops = h_new();

//     for (k=1; k<=Ndata; ++k) {
//         Yk = h_get(Y, Yk_list(k));

//         Ck = D_single_3D(X, Yk, mode, status=status);
//         if(is_void(Ck)) {
//             if (flag_disp_status) trtt_error_display_msg, status=status;
//             return;
//         }

//         h_set, Cops, trtt_get_key_name(k, "C"), Ck;
//     }

//     Ck_list = h_keys(Cops);  
//     h_set, Cops, Ck_list = Ck_list(sort(Ck_list));
    
//     return Cops;
// }
// D_whole_3D = trtt_whole_distance_driven_3D_projector;

/*|------------------------|*/
/*| 3D COMBINED PROJECTORS |*/
/*|------------------------|*/

func trtt_single_3D_combined_projector(X, Yk, mode, type, status=)
/* DOCUMENT trtt_single_3D_combined_projector, X, Yk, mode, type, status=

   Create a 3D+t projector, combining a 3D spatial projector with the model
   given by TYPE (TRTT_SPLINE_DRIVEN or TRTT_DISTANCE_DRIVEN), and a time
   interpolator, according to the geometry of projection MODE
   (TRTT_PARALLEL_BEAM or TRTT_FAN(CONE)_BEAM), to perform Dynamic tomography. X
   is the TOMOGRAPHIC OBJECT and Yk is the DATA object. The projector consists in
   a h_evaluator carrying references to the projector Ck (carrying references to
   X and Yk) and and the time interpolator Bk, associated with the projection
   function. A STATUS can be given. If not it is created in the function and if
   an error occurs, it will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, _trtt_single_2D_combined_projector, h_evaluator.
 */  
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    /* Yk must be a tomdata */
    if (!is_tomdata(Yk)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return;
    }

    if (type == TRTT_SPLINE_DRIVEN) {
        Ak = A_single_3D(X, Yk, mode, status=status);
        if(is_void(Ak)) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    }
    // else if (type == TRTT_DISTANCE_DRIVEN) {
        // Ak = D_single_3D(X, Yk, mode, status=status);
        // if(is_void(Ak)) {
        //     if (flag_disp_status) trtt_error_display_msg, status=status;
        //     return;
        // }
    // }

    Bk = B_single_3D(X, Yk, status=status);
    if(is_void(Bk)) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    p = h_new(Ak = Ak,
              Bk = Bk);

    h_evaluator, p, "_trtt_single_3D_combined_projector";

    return p;
}
C_single_3D = trtt_single_3D_combined_projector;

func _trtt_single_3D_combined_projector(p, x, job)
/* DOCUMENT _trtt_single_3D_combined_projector, p, x, job

   Evaluator (direct and tranpose) of a 3D combined projector. The function
   applies a combination of a spatial 3D projector Ck and a time interpolator
   Bk.

   JOB=0 (direct)

        out = Ck[Bk(x)]

   JOB=1 (transpose)

        out = Bk[Ck(x,1),1]

   SEE ALSO: trtt_single_2D_combined_projector, h_evaluator.
 */ 
{
    Ak = p.Ak;
    Bk = p.Bk;

    if (!job) {
        out = Ak(Bk(x));
    } else if (job==1) {
        out = Bk(Ak(x,1),1);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    return out;
}

func trtt_whole_3D_combined_projector(X, Y, mode, type, status=)
/* DOCUMENT trtt_whole_3D_combined_projector, X, Y, mode, type, status=

   Create a set of 3D combined projectors, from a TOMOGRAPHIC OBJECT X and a set
   of DATA objects Y, according to the geometry of projection MODE. A STATUS can
   be given. If not it is created in the function and if an error occurs, it
   will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, trtt_single_2D_combined_projector.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    // FIXME: Check if Y is a hash table containing TOMDATA objects
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return;
    }
    
    Yk_list = Y.Yk_list;
    ndata = numberof(Yk_list);

    /* Create then hash table which contains the whole projectors */
    Cops = h_new();

    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y, Yk_list(k));

        Ck = C_single_3D(X, Yk, mode, type, status=status);
        if(is_void(Ck)) {
            trtt_error_display_msg, status=status;
            return;
        }

        h_set, Cops, trtt_get_key_name(k, "C"), Ck;
    }

    Ck_list = h_keys(Cops);
    h_set, Cops, Ck_list = Ck_list(sort(Ck_list));
    
    return Cops;
}
C_whole_3D = trtt_whole_3D_combined_projector;
