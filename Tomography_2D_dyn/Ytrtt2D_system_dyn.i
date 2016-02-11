/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_2D_system.i --
 *
 * TRTT 2D simulated system package.
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

require, "random.i";

/* GLOBALS =================================================================== */

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

/* USER FUNCTIONS ============================================================ */

func trtt_2D_create_whole_simu_system(data,
                                      nv,
                                      v_scl,
                                      nx,
                                      ny,
                                      nt,
                                      s_scl,
                                      tcycl,
                                      tacq,
                                      s_deg,
                                      t_deg,
                                      size_footprint,
                                      ndata,
                                      mode,
                                      type,
                                      Rsc,
                                      Rsd,
                                      dates=,
                                      xoff=,
                                      yoff=,
                                      ref_obj=,
                                      voff=,
                                      theta=,
                                      var_noise_uniform=,
                                      photon_flux=,
                                      read_noise=,
                                      add_noise=,
                                      matrix=,
                                      use_sparse_coefs=,
                                      epsilon=,
                                      max_doublings=,
                                      status=)
/* DOCUMENT trtt_2D_create_whole_simu_system, data,
                                              nv,
                                              v_scl,
                                              nx,
                                              ny,
                                              nt,
                                              s_scl,
                                              tcycl,
                                              tacq,
                                              s_deg,
                                              t_deg,
                                              size_footprint,
                                              ndata,
                                              mode,
                                              type,
                                              Rsc,
                                              Rsd,
                                              dates=,
                                              xoff=,
                                              yoff=,
                                              ref_obj=,
                                              voff=,
                                              theta=,
                                              var_noise_uniform=,
                                              photon_flux=, 
                                              read_noise=,
                                              add_noise=,
                                              matrix=,
                                              use_sparse_coefs=,
                                              epsilon=,
                                              max_doublings=
                                              status=

   Create a whole 2D simulation system for tomographic reconstruction with TRTT
   code. The TOMOGRAPHIC OBJECT and the DATA are based on a set of ellipses
   ELLIPSES defining the object of interest (Shepp Logan phantom).
   
   ELLISPES:             the set of ellipses.
   |
   |
   |----- E1
   |       |----- x0     // x-coordinate of the center in meters
   |       |----- y0     // y-coordinate of the center in meters
   |       |----- a      // length of the semimajor axis in meters
   |       |----- b      // length of the semiminor axis in meters
   |       |----- alpha  // angle of rotation of the major axis from
                         // the x-axis in radians
   |       |----- mu     // attenuation (constant on the surface of
                         // the ellipse)
   |
   |
   |----- E2
   |       |----- x0        
   |       |----- y0       
   |       |----- a         
   |       |----- b         
   |       |----- alpha   
   |       |----- mu      
   |
   |
   |----- E3
   :
   :  
   :

   NV:                   detector size.
   V_SCL:                pixels detector sampling (in cm).
   NX:                   x-dimesion of the object of interest.
   NY:                   y-dimesion of the object of interest.
   S_SCL:                spatial sampling step of the object of interest (in cm).
   S_DEG:                spatial spline degree for the model.
   SIZE_FOOTPRINT:       number of points of the footprint profile.
   NDATA:                number of source positions (= number of data objects).
   MODE:                 projection mode (TRTT_PARALLEL_BEAM or TRTT_FAN(CONE)_BEAM).
   TYPE:                 type of projector (TRTT_SPLINE_DIVEN or TRTT_DISTANCE_DRIVEN).
   VOFF:                 voff of the detector (in cm).
   RSC:                  source/center distance (in cm).
   RSD:                  focale (source/detector orthogonal distance) (in cm).
   THETA=:               a set of projection angle.
                         By default NDATA regularly spaced angles
                         are determined on PI range.                       
   PHOTON_FLUX=:         Number of detected photons without absorption (for noise
                         calculation).
   READ_NOISE=:          Pixels detector reading noise (in e-/pixel/image).  
   MATRIX=:              if =1, calculate the projector matrix coefficients.
   EPSILON=:             precision of Romberg integrale
                         (used for analytic sinogram calculation).
   MAX_DOUBLINGS=:       dichotomy depth for Romberg integrale.
   STATUS=:              error controler. If not set, it is created in the function
                         and if an error occurs, it will be directly displayed.

   Such an object definition allows to calculate analytic sinograms and get a
   reference image. Hence from ELLIPSES and all the function arguments, we
   create the reference tomographic object XREF, the tomographic object X used
   for reconstruction which carries all its parameters, and the set of data objects
   Y which carries the sinogram data and the calibration parameters.
   _
    |
    |_ X             => TOMOBJ
    |
    |_ XREF          => TOMOBJ
    |
    |_ Y
       |
       |_ Y01        => TOMDATA
       |
       |_ Y02        => TOMDATA
       |
       |_ Yk         => TOMDATA
       |
       :
       :
       |_ Yndata     => TOMDATA
       |
       |_ Yk_list = ["Y01","Y02",...,"Yk",...,"Yndata"]
   
   From these features, the set of projectors COPS created from X and each Yk.
   _
    |
    |_ COPS
       |
       |_ C01        => PROJECTOR
       |
       |_ C02        => PROJECTOR
       |
       |_ Ck         => PROJECTOR
       |
       :
       :
       |_ Cndata     => PROJECTOR
       |
       |_ Ck_list = ["C01","C02",...,"Ck",...,"Cndata"]

   If MATRIX=1, the coefficients of the whole projection matrix are calculated
   and stored as a sparse matrix.

   A hashtab HTOMO with all these entries is returned.

   SEE ALSO: trtt_tomobj_info, trtt_tomdata_info, trtt_2D_projectors_info,
             trtt_error_info.
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

    /* Create the set of TOMDATA */
    Y = trtt_2D_create_simu_tomdata_set(data, nv, v_scl, ndata, nt, tcycl, tacq, mode, dates=dates, xoff=xoff, yoff=yoff, voff=voff, theta=theta, Rsc=Rsc, Rsd=Rsd, epsilon=epsilon, max_doublings=max_doublings, status=status);
    if (is_void(Y)) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* FIXME: Data and noise covariance determination for weighting the data
     * residuals */
    
    /* Deal with uniform noise if VAR_NOISE_UNIFORM specified */
    if (!is_void(var_noise_uniform)) {
        if (trtt_2D_deal_uniform_noise(Y, var_noise_uniform, add_noise=add_noise, status=status) == TRTT_FAILURE) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    }
    /* Deal with non uniform noise if PHOTON_FLUX specified */
    else if (!is_void(photon_flux)) {
        if (trtt_2D_deal_non_uniform_noise(Y, photon_flux, read_noise=read_noise, add_noise=add_noise, status=status) == TRTT_FAILURE) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    }
    /* Deal with no noise specified */
    else {
        if (trtt_2D_deal_no_noise(Y, status=status) == TRTT_FAILURE) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    }

    /* Create the reference TOMOBJ */
    if (!is_void(voxels_ref)) {
        Xref = trtt_2D_create_tomobj(nx, ny, nt, s_scl, tcycl, s_deg, t_deg, size_footprint, voxels_ref=ref_obj, status=status);
        if (is_void(Xref)) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    }

    /* Create the reconstruction TOMOBJ */
    X = trtt_2D_create_tomobj(nx, ny, nt, s_scl, tcycl, s_deg, t_deg, size_footprint, status=status);
    if (is_void(X)) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    /* Create the projectors */
    Cops = trtt_2D_create_projectors(X, Y, mode, type, use_sparse_coefs=use_sparse_coefs, matrix=matrix, status=status);
    
    /* Create the whole system hash tab */
    return h_new(X=X, Y=Y, Xref=Xref, Cops=Cops);
}
trtt_2D_create = trtt_2D_create_whole_simu_system;

func trtt_2D_load_simu_system(filename, matrix=)
/* DOCUMENT trtt_2D_load_simu_system, filename, matrix=

   Load a whole 2D system for tomographic reconstruction with TRTT code from a
   YHD file FILENAME. All the entries, the tomographic object X, the reference
   object Xref, the set of data Y and the projectors COPS are recovered and
   returned in a hashtab. If MATRIX=1, and if it is not done yet, the sparse
   projection matrix is calculated.

   SEE ALSO:  trtt_2D_create_whole_simu_system, trtt_2D_save_simu_system,
              yhd_restore.
 */ 
{
    Htomo = yhd_restore(filename);

    X = Htomo.X;
    Xref = Htomo.Xref;
    Y = Htomo.Y;
    Cops = Htomo.Cops;
    eval = Cops.eval;
    evalAk = Cops.evalAk;
    evalBk = Cops.evalBk;    

    /* Re create xforms for TOMOBJ */
    obj_xi_coefs = X.xi_coefs;
    obj_xi = trtt_xform3d_new();
    trtt_xform3d_set_coefficients, obj_xi, obj_xi_coefs;
    h_set, X, xi = obj_xi;

    /* Re create interp_coefs for TOMOBJ */
    nx = X.nx;
    ny = X.ny;
    nz = X.nz;
    nt = X.nt;
    s_deg = X.s_deg;
    t_deg = X.t_deg;
    /* Create interp coefs */
    wx = spl_interp_coefs_new_spline(indgen(nx)-1.0, nx, s_deg);
    wy = spl_interp_coefs_new_spline(indgen(ny)-1.0, ny, s_deg);
    wz = spl_interp_coefs_new_spline(indgen(nz)-1.0, nz, s_deg);
    wt = spl_interp_coefs_new_spline(indgen(nt)-1.0, nt, t_deg);
    /* Set interp coefs */
    h_set, X, wx=wx;
    h_set, X, wy=wy;
    h_set, X, wz=wz;
    h_set, X, wt=wt;
    if (!is_void(Xref)) {
        h_set, Xref, wx=wx;
        h_set, Xref, wy=wy;
        h_set, Xref, wz=wz;
        h_set, Xref, wt=wt;
    }
    
    /* Re create xforms for TOMDATA */
    Yk_list = Y.Yk_list;
    Ck_list = Cops.Ck_list;
    ndata = numberof(Yk_list);
    nv = h_get(Y,Yk_list(1)).nv;
    
    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y,Yk_list(k));
        Ck = h_get(Cops,Ck_list(k));
        Ak = Ck.Ak;
        Bk = Ck.Bk;
        /* Restore evaluators */
        h_evaluator, Ck, eval;
        h_evaluator, Ak, evalAk;
        h_evaluator, Bk, evalBk;
        /* Get coefs */
        data_xi_coefs = Yk.xi_coefs;
        xotilt_coefs = Yk.xotilt_coefs;
        xos_coefs = Yk.xos_coefs;
        xsd_coefs = Yk.xsd_coefs;
        /* Create xforms */
        data_xi = trtt_xform3d_new();
        xotilt = trtt_xform3d_new();
        xos = trtt_xform3d_new();
        xsd = trtt_xform3d_new();
        /* Set coefs */
        trtt_xform3d_set_coefficients, data_xi, data_xi_coefs;
        trtt_xform3d_set_coefficients, xotilt, xotilt_coefs;
        trtt_xform3d_set_coefficients, xos, xos_coefs;
        trtt_xform3d_set_coefficients, xsd, xsd_coefs;
        /* Set xforms */
        h_set, Yk, xi = data_xi;
        h_set, Yk, xotilt = xotilt;
        h_set, Yk, xos = xos;
        h_set, Yk, xsd = xsd;

        /* Re hang projectors parameters */
        h_set, Ak, X=X;
        h_set, Ak, Yk=Yk;

        /* Create and set time interpolation coefficients */
        t_index = Yk.t_index;
        wtBk = spl_interp_coefs_new_spline(t_index, nt, t_deg);
        h_set, Bk, wt=wtBk;
    }
    /* Set projector matrix if it do not exist */
    if (matrix == TRTT_TRUE) {
        _trtt_2D_sparse_matrix, Cops;
    }

    return Htomo;
}
trtt_2D_load = trtt_2D_load_simu_system;

func trtt_2D_save_simu_system(Htomo, filename, overwrite=)
/* DOCUMENT trtt_2D_save_simu_system, Htomo, filename, overwrite=

   Save a whole 2D system for tomographic reconstruction with TRTT code in a YHD
   file FILENAME. All entries are re-arranged to store only useful
   parameters. When such a file is loaded with the function
   trtt_2D_load_simu_system(), all entries are recovered in its original
   state. If keyword OVERWRITE is true and file FILENAME already exists, the new
   file will (silently) overwrite the old one; othwerwise, file FILENAME must
   not already exist (default behaviour).

   SEE ALSO:  trtt_2D_create_whole_simu_system, trtt_2D_load_simu_system,
              yhd_save.
 */ 
{
    Htomo_temp = h_copy(Htomo);
    X = Htomo_temp.X;
    Xref = Htomo_temp.Xref
    Y = Htomo_temp.Y;
    Cops = Htomo_temp.Cops;

    Yk_list = Y.Yk_list;
    Ck_list = Cops.Ck_list;
    ndata = numberof(Yk_list);

    /* Remove unsavable keys */
    h_delete, X, "wx";
    h_delete, X, "wy";
    h_delete, X, "wz";
    h_delete, X, "wt";
    h_delete, X, "xi";
    if (!is_void(Xref)) {
        h_delete, Xref, "wx";
        h_delete, Xref, "wy";
        h_delete, Xref, "wz";
        h_delete, Xref, "wt";
        h_delete, Xref, "xi";
    }

    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y,Yk_list(k));
        Ck = h_get(Cops,Ck_list(k));
        Ak = Ck.Ak;
        Bk = Ck.Bk;
        h_delete, Yk, "xi";
        h_delete, Yk, "xotilt";
        h_delete, Yk, "xos";
        h_delete, Yk, "xsd";
        h_delete, Ak, "Yk";
        h_delete, Ak, "X";
        if (!is_void(Ak.MCk))
            h_delete, Ak, "MCk";
        h_delete, Bk, "wt";
    }
    
    yhd_save, filename, Htomo_temp, overwrite=overwrite;
}
trtt_2D_save = trtt_2D_save_simu_system;

func trtt_2D_create_simu_tomdata_set(data, nv, v_scl, ndata, nt, tcycl, tacq, mode, dates=, xoff=, yoff=, voff=, theta=, Rsc=, Rsd=, epsilon=, max_doublings=, status=)
/* DOCUMENT trtt_2D_create_simu_tomdata_set, data, nv, v_scl, ndata, nt, tcycl, tacq, mode, dates=, xoff=, yoff=, voff=, theta=, Rsc=, Rsd=, epsilon=, max_doublings=, status=

   Create a set of DATA objects Y of a whole 2D simulation system for
   tomographic reconstruction with TRTT code, and calculate the analytic
   sinograms. A STATUS can be given. If not it is created in the function and if
   an error occurs, it will be directly displayed.
   
   SEE ALSO: trtt_2D_create_whole_simu_system, trtt_tomdata_info,
             trtt_ellipses_fan_beam_analytic_Radon_transform,
             trtt_ellipses_parallel_beam_analytic_Radon_transform.                        
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
    
    if (mode != TRTT_PARALLEL_BEAM & mode != TRTT_FAN_BEAM & mode != TRTT_CONE_BEAM) {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    if ((mode==TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM)
        & (is_void(Rsc) || is_void(Rsd))) {
        trtt_error_set, status, TRTT_ERR_MISSING_PARAMETER,
            msg="You must specify the source/center and source/detector distances in fan beam mode.",
            disp=flag_disp_status;
        return;
    }

    /* "Invisible" parameters */
    nu = 1;
    u_scl = 1.0;

    /* Determine time indexes for projections */
    if (is_void(dates)) {
        t = indgen(0:ndata-1)*(tacq/(ndata-1.));
        t_index = (t%tcycl)*(nt/tcycl);
    } else {
        if (numberof(dates)!=ndata) {
            trtt_error_set, status, TRTT_ERR_BAD_NDIMS,
                msg="There must be as DATES as the number of projections NDATA.",
                disp=flag_disp_status;
        } else {
            tacq=dates(0);
            local t; eq_nocopy, t, dates;
            t_index = (t%tcycl)*(nt/tcycl);
        }
    }

    /* Projection angles */
    // FIXME: the vector of angles THETA can be directly specified in function
    // parameters. Else it is calculated from NDATA.
    if (is_void(theta)) {
        theta_step = pi/ndata;
        theta_range = theta_step*(ndata-1);
        theta=trtt_span(0.0, theta_range, ndata);
    }
    
    if (!is_void(voff) & (is_array(voff) || is_scalar(voff))) {
        if (numberof(voff) != ndata & numberof(voff) != 1)
            trtt_error_set, status, TRTT_ERR_BAD_NDIMS,
                msg="Offset parameter VOFF has bad number of elements.",
                disp=flag_disp_status;
    } 
    

    /* Create Data hash tab */
    Y = h_new();
    // FIXME: Each single TOMDATA is created ans inserted in Data hash tab
    for (k=1; k<=ndata; ++k) {
        /* Create TOMDATA */
        Yk = trtt_tomdata_new(data(,k), u_scl, v_scl, t_index(k), status=status);
        /* Extract xforms */
        xos = Yk.xos;
        xsd = Yk.xsd;
        xotilt = Yk.xotilt;
        
        /* Quaternion : rotation of theta around z axis */
        th_k = 0.5*(pi-theta(k));
        qw = cos(th_k);
        qx = qy = 0.0;
        qz = sin(th_k);
        trtt_xform3d_rotate, xos, qw, qx, qy, qz, xos;
        
        trtt_xform3d_move, xos, Rsc, 0.0, 0.0, xos;
        trtt_xform3d_move, xsd, -1.0*Rsd, 0.0, 0.0, xsd;

        if (!is_void(voff)) {
            if (!is_scalar(voff))
                trtt_xform3d_move, xsd, 0.0, voff(k), 0.0, xsd;
            else
                trtt_xform3d_move, xsd, 0.0, voff, 0.0, xsd;
        }

        if (!is_void(xoff))
            trtt_xform3d_move, xotilt, xoff, 0.0, 0.0, xotilt;

        if (!is_void(yoff))
            trtt_xform3d_move, xotilt, 0.0, yoff, 0.0, xotilt;
        
        /* Save xform coefficients */
        h_set, Yk, xi_coefs = Yk.xi.a;
        h_set, Yk, xotilt_coefs = Yk.xotilt.a;
        h_set, Yk, xos_coefs = Yk.xos.a;
        h_set, Yk, xsd_coefs = Yk.xsd.a;

        
        h_set, Yk, theta=theta(k);
        /* Set the current TOMDATA in Data hash tab */
        h_set, Y, trtt_get_key_name(k, "Y"), Yk;
    }

    Yk_list = h_keys(Y);
    h_set, Y, Yk_list = Yk_list(sort(Yk_list));

    return Y;
}

func trtt_2D_create_tomobj(nx, ny, nt, s_scl, tcycl, s_deg, t_deg, size_footprint, voxels_ref=, status=)
/* DOCUMENT trtt_2D_create_tomobj, nx, ny, nt, s_scl, tcycl, s_deg, t_deg,
   size_footprint, voxels_ref=, status=

   Create a TOMOGRAPHIC OBJECT X of a whole 2D simulation system for tomographic
   reconstruction with TRTT code. If VOXELS_REF is given, then it is a reference
   object and the true voxels image is inserted. A STATUS can be given. If not
   it is created in the function and if an error occurs, it will be directly
   displayed.
   
   SEE ALSO: trtt_2D_create_whole_simu_system, trtt_tomobj_info,
             trtt_create_ellipses_phantom. 
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

    /* "Invisible" parameters */
    nz = 1;
    /* temporal step */
    t_scl = tcycl/nt;

    /* Create tomobj */
    X = trtt_tomobj_new(array(double,nx,ny,nt), s_scl, t_scl, s_deg, t_deg, size_footprint, status=status);
    /* Save xform coefficients */
    h_set, X, xi_coefs = X.xi.a;

    if (!is_void(voxels_ref)) {
        /* Create voxels image */
        trtt_tomobj_set_voxels, X, voxels_ref;
    }

    return X;
}

func trtt_2D_create_projectors(X, Y, mode, type, use_sparse_coefs=, matrix=, status=)
/* DOCUMENT trtt_2D_create_projectors, X, Y, mode, type, use_sparse_coefs=, matrix=, status=

   Create a set of 2D projectors COPS of a whole 2D simulation system for
   tomographic reconstruction with TRTT code. Two types of projectors are
   available: the SEPARABLE B-SPLINE projector and the DISTANCE DRIVEN projector
   ; and for two different geometries of projection: PARALLEL BEAM and FAN/CONE
   BEAM. If MATRIX=1, the sparse projection matrix is calculated. A STATUS can
   be given. If not it is created in the function and if an error occurs, it
   will be directly displayed.
   
   SEE ALSO: trtt_2D_create_whole_simu_system, trtt_2D_projectors_info. 
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

    if (mode != TRTT_PARALLEL_BEAM & mode != TRTT_FAN_BEAM & mode != TRTT_CONE_BEAM) {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    if (type != TRTT_SPLINE_DRIVEN & type != TRTT_DISTANCE_DRIVEN) {
        trtt_error_set, status, TRTT_ERR_BAD_PROJECTOR, disp=flag_disp_status;
        return;
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    /* FIXME: Check if Y is a hash table containing TOMDATA objects */
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return;
    }

    Cops = C_whole_2D(X, Y, mode, type, status=status);

    h_set, Cops, eval="_trtt_single_2D_combined_projector";
    h_set, Cops, evalAk = "_trtt_single_2D_projector";
    h_set, Cops, evalBk = "_trtt_single_2D_time_interpolator";
    h_set, Cops, type = type;

    if (matrix==TRTT_TRUE) {
        _trtt_2D_sparse_matrix, Cops, dyn=1n;
    }

    if (use_sparse_coefs==TRTT_TRUE) {
        _trtt_2D_sparse_coefs, Cops, dyn=1n;
    }

    return Cops;
}

func trtt_2D_deal_no_noise(Y, status=)
/* DOCUMENT trtt_2D_deal_no_noise, Y, status=

   Calculate weights of inverse data covariance matrix, considering a the global
   data variance VAR_DATA, from a set of DATA objects Y.  The function stores in
   each Yk the array of weights v(j,k) for noise covariance matrix
   calculation. As variables are independent, this matrix will be diagonal.

   A STATUS can be given. If not it is created in the function and if an error
   occurs, it will be directly displayed.
 */  
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) { /* Check if there is already an error */
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init(void);
    }

    /* FIXME: Check if Y is a hash table containing TOMDATA objects */
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    Yk_list = Y.Yk_list;
    ndata = numberof(Yk_list);

    /* Get noise-free sinogram */
    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y,Yk_list(k));
        nv = Yk.nv;
        /* Covariance weights */
        h_set, Yk, weight = array(1.0,nv);
        h_set, Yk, mean_var_noise = 0.0;
    }
    
    return TRTT_SUCCESS;
}

func trtt_2D_deal_uniform_noise(Y, var_noise, add_noise=, status=)
/* DOCUMENT trtt_2D_deal_uniform_noise, Y, var_noise, add_noise=, status=

   Calculate weights of inverse noise covariance matrix, considering a
   stationary gaussian noise, with variance VAR_NOISE, from a set of DATA
   objects Y.  The function stores in each Yk the array of weights v(j,k) for
   noise covariance matrix calculation. As variables are independent, this
   matrix will be diagonal.

   If ADD_NOISE is set to 1, a noise, with the calculated variance, is added to
   the noiseless data in each Yk (this case is for numerically simulated data).

   A STATUS can be given. If not it is created in the function and if an error
   occurs, it will be directly displayed.
 */  
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) { /* Check if there is already an error */
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init(void);
    }

    /* FIXME: Check if Y is a hash table containing TOMDATA objects */
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    Yk_list = Y.Yk_list;
    ndata = numberof(Yk_list);

    /* Get noise-free sinogram */
    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y,Yk_list(k));
        nv = Yk.nv;
        /* Covariance weights */
        inv_var_noise = 1.0/var_noise;
        h_set, Yk, weight = array(inv_var_noise,nv);
        h_set, Yk, mean_var_noise = var_noise;
        /* Noisy projection */
        if(add_noise) {
            pixels = Yk.pixels;
            // noisy_pixels = max(0.0,pixels + random_n(nv)*sqrt(var_noise));
            noisy_pixels = pixels + random_n(nv)*sqrt(var_noise);
            trtt_tomdata_set_pixels, Yk, noisy_pixels;
        }
    }
    
    return TRTT_SUCCESS;
}


func trtt_2D_deal_non_uniform_noise(Y, photon_flux, read_noise=, add_noise=, status=)
/* DOCUMENT trtt_2D_deal_non_uniform_noise, Y, photon_flux, read_noise=, add_noise=, status=

   Calculate weights of inverse noise covariance matrix, considering a
   non-stationary gaussian noise, from a set of DATA objects Y. PHOTON_FLUX is
   the number of detected photons without absorption and READ_NOISE is the
   pixels detector reading noise in e-/pixel/image. These parameters are used to
   calculate the noise variance which follows the law :
   
                         z(j,k) + var_CCD
   var_noise = v(j,k) =  ---------------- with z(j,k) = lambda_0*exp[-q(j,k)]
                             (z(j,k))²
   
   Where q(j,k) is the projection of the attenuation on pixel detector j,
   integrated on its response, at projection angle k. LAMBDA_0 is the photon
   flux and var_CCD is the variance of the reading noise.

   The noise is supposed to be the sum of a Poisson noise (attenuation of a
   photon flux : Beer-Lambert's law) and a gaussian reading noise on the
   detector. We make the assumption of low absorption, from which we deduce that
   the measures y(j,k) are independent gaussian random variables (non stationary
   gaussian noise).
   
              /                 \
   y(j,k) = N | q(j,k) , v(j,k) |
              \                 /

   See http://mitiv.univ-lyon1.fr/index.php/Tomographie/NoiseModels
   for the detailed noise study.
                     
   The function stores in each Yk the array of weights v(j,k) for noise
   covariance matrix calculation. As variables are independent, this matrix will
   be diagonal.

   If ADD_NOISE is set to 1, a noise, with the calculated variance, is added to
   the noiseless data in each Yk (this case is for numerically simulated data).

   A STATUS can be given. If not it is created in the function and if an error
   occurs, it will be directly displayed.
 */  
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) { /* Check if there is already an error */
            return TRTT_FAILURE;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init(void);
    }

    /* FIXME: Check if Y is a hash table containing TOMDATA objects */
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return TRTT_FAILURE;
    }

    Yk_list = Y.Yk_list;
    ndata = numberof(Yk_list);

    if (is_void(read_noise)) read_noise = 0.0;

    /* Get noise-free sinogram */
    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y,Yk_list(k));
        nv = Yk.nv;
        /* Noise variance */
        if (photon_flux != 0.0) {
            pixels = Yk.pixels;
            z_jk = photon_flux * exp(-1.0*pixels);
            z_jk_sqr = z_jk * z_jk;
            var_noise = (z_jk + read_noise^2.0)/z_jk_sqr;
            mean_var_noise = (1./nv)*sum(var_noise(*));
        } else {
            var_noise = 1.0;
            mean_var_noise = 0.0;
        }
        h_set, Yk, weight = 1.0/(var_noise(*));
        h_set, Yk, mean_var_noise = mean_var_noise;
        
        /* Noisy projection */
        if(add_noise) {
            noisy_pixels = max(0.0,pixels + random_n(nv)*sqrt(var_noise));
            trtt_tomdata_set_pixels, Yk, noisy_pixels;
        }
    }
    
    return TRTT_SUCCESS;
}

