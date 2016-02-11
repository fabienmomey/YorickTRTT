/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_2D_projectors_comparator.i --
 *
 * TRTT 2D Radon projectors for verifying the used operators.
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

/* GLOBALS =================================================================== */

/* INTERNAL LOW LEVEL FUNCTIONS ============================================== */

/*
 * This function performs the same integration algorithm as Munro's one for
 * the INTEG function in Yorick. Y is the array of M points to integrate from
 * its initial point to NP regularly spaced bounds, starting from XP_I_START
 * (expressed in index units of Y), with an increment INC between points. YCUM
 * is a M-points array giving the partial sums of Y, considering a sampling
 * STEP between points.
 */
func _trtt_integ_munro(yp, y, ycum, xp_i_start, inc, np, m, step)
{
  /* Compute interpolated integrals */
  for (i=0; i<np; ++i) {
    xp_i = xp_i_start + i*inc;
    ip = long(floor(xp_i));
    if (ip>=0 & ip<m-1) {
      dx = xp_i-ip;
      c1 = 0.5*dx*dx;
      c0 = dx-c1;
    }
    
    if (ip<0) {
        yp(i+1) = 0.0; 
    } else if (ip<m-1) {
        yp(i+1) = step*(c0*y(ip+1)+c1*y(ip+2))+ycum(ip+1);
    } else {
        yp(i+1) = ycum(m);
    } 
  }
}

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

func _trtt_2D_projection_O_to_D_spline_driven (xn, yn, xos, xsd, mode)
{
    /* O to S transform */
    xs = xos(1)*xn + xos(2)*yn(-,) + xos(4);
    ys = xos(5)*xn + xos(6)*yn(-,) + xos(8);

    /* S to D projection */
    if (anyof(xs==0.0) & mode == TRTT_FAN_BEAM)
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";

    a1 = xsd(1);
    if (a1 == 0)
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";
    
    a1inv = 1.0/a1;
    a2 = xsd(2);
    a4 = xsd(4);
    a5 = xsd(5);
    a6 = xsd(6);
    a8 = xsd(8);

    b1 = a6 + a1inv*a5*a2;
    b3 = a8 + a1inv*a5*a4;

    if (mode == TRTT_PARALLEL_BEAM) {
        return (b1*ys+b3);
    } else {
        tan_gamma = ys/xs;
        xeq = -(1.0/(a1+a2*tan_gamma))*a4; 
        return (b1*tan_gamma*xeq + b3);
    }
}

func _trtt_y2D_projection_O_to_S(coord, xos)
{
    outcoord = array(double,2);
    x = coord(1);
    y = coord(2);
    /* O to S transform */
    outcoord(1) = xos(1)*x + xos(2)*y + xos(4);
    outcoord(2) = xos(5)*x + xos(6)*y + xos(8);

    return outcoord;
}

func _trtt_y2D_projection_S_to_D(coord, xsd, mode)
{
    outcoord = array(double,2);
    xs = coord(1);
    ys = coord(2);
    
    /* S to D projection */
    if (anyof(xs==0.0) & mode == TRTT_FAN_BEAM)
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";

    a1 = xsd(1);
    if (a1 == 0)
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";
    
    a1inv = 1.0/a1;
    a2 = xsd(2);
    a4 = xsd(4);
    a5 = xsd(5);
    a6 = xsd(6);
    a8 = xsd(8);

    b1 = a6 + a1inv*a5*a2;
    b3 = a8 + a1inv*a5*a4;

    if (mode == TRTT_PARALLEL_BEAM) {
        outcoord(1) = a1inv*(a2*ys + a4);
        outcoord(2) = b1*ys + b3;
        return outcoord;
    } else {
        tan_gamma = ys/xs;
        xeq = -(1.0/(a1+a2*tan_gamma))*a4;
        outcoord(1) = xeq;
        outcoord(2) = b1*tan_gamma*xeq + b3;
        return outcoord;
    }
}

func _trtt_y2D_burn_footprint(pix, vox, coord, v_scl, nv, s_scl, s_deg, rho_0, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint)
{
  /* Normalized b-spline step to pixel detector sampling  */
  h_norm = (1.0/v_scl)*(rho_0/s_scl);
  scl_norm = 1/(h_norm*step_footprint);
  /* Retrieve the half support of the footprint on detector */
  half_support = 0.5 * (s_deg+1) * h_norm;
  /* Get the position of impact point of the blob center */
  vc =  coord(2);
  /* Get the bounds positions of the footprint */
  v_left = vc-half_support;
  v_right = vc+half_support;
  // FIXME: we are in bound detector pixel coordinates
  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv

  /* Get the indexes of impacted pixels detector */
  d_left = floor(v_left);
  d_right = floor(v_right);
  d_left = min(max(d_left,0),nv-1);
  d_right = max(min(d_right,nv-1),0);
  
  if (d_left <= d_right) {
    /* Get the position of left-hand side of the footprint in pixel detector
     * space */
      coord_footprint_0 = vc + coord_footprint(1)*h_norm;
      coord_footprint_start = (double(d_left)-coord_footprint_0)*scl_norm;
      // FIXME: footprint samples coordinates
      // FIXME: 0 ; 1 ; 2 ; ... ; size_footprint-1
    
      npix = long(d_right-d_left+1); // FIXME: number of impacted detector pixels

      // FIXME: bound effects are handled by INTEG function
      
      pix_integ = array(double, npix+1);
    
      _trtt_integ_munro, pix_integ, footprint,                           \
          cum_footprint, coord_footprint_start,                         \
          scl_norm, npix+1, size_footprint, step_footprint;

      for (i=1; i<=npix; ++i) {
          pix(long(d_left)+i) += rho_0*vox*(pix_integ(i+1)-pix_integ(i));
      }    
  }
}

func _trtt_y2D_burn_footprint_transp(vox, pix, coord, v_scl, nv, nx, ny, id_x, id_y, s_scl, s_deg, rho_0, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint)
{
  /* Normalized b-spline step to pixel detector sampling  */
  h_norm = (1.0/v_scl)*(rho_0/s_scl);
  scl_norm = 1/(h_norm*step_footprint);
  /* Retrieve the half support of the footprint on detector */
  half_support = 0.5 * (s_deg+1) * h_norm;
  /* Get the position of impact point of the blob center */
  vc =  coord(2);
  /* Get the bounds positions of the footprint */
  v_left = vc-half_support;
  v_right = vc+half_support;
  // FIXME: we are in bound detector pixel coordinates
  // FIXME: 0 ; 1 ; 2 ; ... ; nv-1 ; nv

  /* Get the indexes of impacted pixels detector */
  d_left = floor(v_left);
  d_right = floor(v_right);
  d_left = min(max(d_left,0),nv-1);
  d_right = max(min(d_right,nv-1),0);
  
  if (d_left <= d_right) {
    /* Get the position of left-hand side of the footprint in pixel detector
     * space */
      coord_footprint_0 = vc + coord_footprint(1)*h_norm;
      coord_footprint_start = (double(d_left)-coord_footprint_0)*scl_norm;
      // FIXME: footprint samples coordinates
      // FIXME: 0 ; 1 ; 2 ; ... ; size_footprint-1
    
      npix = long(d_right-d_left+1); // FIXME: number of impacted detector pixels

    // FIXME: bound effects are handled by INTEG function
      
      pix_integ = array(double, npix+1);
    
      _trtt_integ_munro, pix_integ, footprint,                          \
          cum_footprint, coord_footprint_start,                         \
          scl_norm, npix+1, size_footprint, step_footprint;

      for (i=1; i<=npix; ++i) {
          vox(id_x,id_y) += rho_0 * pix(long(d_left)+i) * (pix_integ(i+1)-pix_integ(i));
      }    
  }
}

func trtt_yPB2D(out, in, nx, ny, s_scl, rho_0, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos, xsd, job)
{
    /*  Calculated parameters */
    coord = array(double,2);
    
    for (j=0; j<ny; ++j) {
        for (i=0; i<nx; ++i) {
            coord(1) = i;
            coord(2) = j;
            
            /* Transform Object -> Source */
            coord = _trtt_y2D_projection_O_to_S(coord, xos);
            /* Projection on detector */
            coord = _trtt_y2D_projection_S_to_D(coord, xsd, TRTT_PARALLEL_BEAM);
            
            if (job == 0) {
                /* DIRECT OPERATOR */
                /* Get voxel value */
                in_ij = in(i+1,j+1);
                _trtt_y2D_burn_footprint, out, in_ij, coord,             \
                    v_scl, nv, s_scl,                                   \
                    s_deg, rho_0,                               \
                    footprint, cum_footprint, size_footprint,           \
                    step_footprint, coord_footprint;
            } else {
                /* TRANSPOSE OPERATOR */
                _trtt_y2D_burn_footprint_transp, out, in, coord, v_scl, nv, \
                    nx, ny, i+1, j+1, s_scl,                                \
                    s_deg, rho_0,                               \
                    footprint, cum_footprint, size_footprint,           \
                    step_footprint, coord_footprint;
            }
        }
    }
}


func trtt_yFB2D(out, x, nx, ny, s_scl, rho_0, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, job)
{
    /*  Calculated parameters */
    coord = array(double,2);
    rho_0_old = rho_0;
    focale = xsd_coefs(4);
    
    for (j=0; j<ny; ++j) {
        for (i=0; i<nx; ++i) {
            coord(1) = i;
            coord(2) = j;
            
            /* Transform Object -> Source */
            coord = _trtt_y2D_projection_O_to_S(coord, xos);
            /* Projection on detector */

            prev_x = coord(1);
            tan_gamma = coord(2)/prev_x;
            
            coord = _trtt_y2D_projection_S_to_D(coord, xsd, TRTT_PARALLEL_BEAM);

            
            if (prev_x != 0) {
                rho_0 = rho_0_old * (focale/prev_x) * sqrt(1.0+tan_gamma*tan_gamma);
            } else {
                fprintf(stderr, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n");
                return;
            }
      
            /* Burn footprint on Detector */
            
            if (job == 0) {
                /* DIRECT OPERATOR */
                /* Get voxel value */
                in_ij = in(i+1,j+1);
                _trtt_y2D_burn_footprint, out, in_ij, coord,             \
                    v_scl, nv, s_scl,                                   \
                    s_deg, rho_0,                               \
                    footprint, cum_footprint, size_footprint,           \
                    step_footprint, coord_footprint;
            } else {
                /* TRANSPOSE OPERATOR */
                _trtt_y2D_burn_footprint_transp, out, in, coord, v_scl, nv, \
                    nx, ny, i+1, j+1, s_scl,                                \
                    s_deg, rho_0,                               \
                    footprint, cum_footprint, size_footprint,           \
                    step_footprint, coord_footprint;
            }
        }
    }
}

/* USER FUNCTIONS ============================================================ */

/*|---------------|*/
/*| 2D PROJECTORS |*/ //FIXME: Same as C-implemented
/*|---------------|*/

func trtt_ysingle_2D_projector(X, Yk, mode, status=)
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
    if (!(t_index>=0 & t_index<=(nt-1))) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
        return;
    }

    p = h_new(X=X,
              Yk=Yk,
              xos_coefs = xos_coefs,
              xsd_coefs = xsd_coefs);

    s_scl = X.s_scl;
    if (mode == TRTT_PARALLEL_BEAM) {
        a5a1 = xsd_coefs(5)/xsd_coefs(1);
        rho_0 = s_scl^2.0 * sqrt(1+a5a1*a5a1);
        h_set, p, rho_0 = rho_0;
        h_set, p, projector = symlink_to_name("trtt_yPB2D")
    }
    else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
        rho_0 = s_scl^2.0;
        h_set, p, rho_0 = rho_0;
        h_set, p, projector = symlink_to_name("trtt_yFB2D")
    }
    else {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    h_evaluator, p, "_trtt_ysingle_2D_projector";

    return p;
}
E_single_2D = trtt_ysingle_2D_projector;

func _trtt_ysingle_2D_projector(p, x, job)
/* DOCUMENT
 */ 
{
    nx = p.X.nx;
    ny = p.X.ny;
    s_scl = p.X.s_scl;
    s_deg = p.X.s_deg;
    rho_0 = p.rho_0;
    local footprint; eq_nocopy, footprint, p.X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, p.X.cum_footprint;
    size_footprint = p.X.size_footprint;
    step_footprint = p.X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, p.X.coord_footprint;
    nv = p.Yk.nv;
    v_scl = p.Yk.v_scl;
    local xos_coefs; eq_nocopy, xos_coefs, p.xos_coefs;
    local xsd_coefs; eq_nocopy, xsd_coefs, p.xsd_coefs;
    __projector = p.projector;

    if (!job) {
        job=0n;
        out = array(double, nv);
    } else if (job==1) {
        out = array(double, nx, ny, 1); // FIXME: nt dimension equals
                                        // 1 but must be specified for
                                        // compatibility with time
                                        // interpolator
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    if(__projector(out, x, nx, ny, s_scl, rho_0, s_deg,
                   footprint, cum_footprint,
                   size_footprint, step_footprint, coord_footprint,
                   nv, v_scl,
                   xos_coefs, xsd_coefs, job) == TRTT_FAILURE) {
        trtt_error_display_msg, code=TRTT_ERR_PROJECTOR_FAILURE;
        return;
    } else {
        return out;
    }             
}

func trtt_ywhole_2D_projector(X, Y, mode, status=)
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

        Ck = E_single_2D(X, Yk, mode, status=status);
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
E_whole_2D = trtt_ywhole_2D_projector;

/*|-------------------|*/
/*| 2D FFT PROJECTORS |*/ //FIXME: Same as Yorick-implemented (using FFT)
/*|-------------------|*/

func trtt_single_spline_driven_2D_projector(X, Yk, mode, status=)
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
    obj_xi = X.xi;
    nt = X.nt;
    s_deg = X.s_deg;
    Np = X.size_footprint;
    /* - Yk */
    nv = Yk.nv;
    v_scl = Yk.v_scl;
    xotilt = Yk.xotilt;
    xos = Yk.xos;
    xsd = Yk.xsd;
    t_index = Yk.t_index;
    /* xform combinations */
    xos_plus = trtt_xform3d_combine(xotilt, obj_xi);
    trtt_xform3d_combine, xos_plus, xos, xos_plus;
    xos_coefs = xos_plus.a;
    xsd_coefs = xsd.a;

    /* Check t_index */
    if (!(t_index>=0 & t_index<=(nt-1))) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
        return;
    }

    /* Calculate the support of the profile */
    support = double(s_deg+1)*sqrt(2*s_scl^2) + v_scl;
    
    /* Profile coordinates */
    step_z = (1.0/(Np-1))*support;
    z = indgen(0:Np-1) - (Np/2) + !(Np%2);
    z *= step_z;
    zmin = min(z);
    zmax = max(z);
    
    /* FFT frequencies coordinates */
    fe = 1.0/step_z;
    w = fft_indgen(Np);
    step_w = fe/Np;
    w *= step_w;

    /* Create the Detector Door for the parallel beam mode */
    if (mode == TRTT_PARALLEL_BEAM) {
        fft_D = v_scl*sinc(w*v_scl);
    }

    p = h_new(X=X,
              Yk=Yk,
              xos_coefs = xos_coefs,
              xsd_coefs = xsd_coefs,
              fft_D = fft_D,
              w = w,
              step_w = step_w,
              z = z,
              step_z = step_z,
              zmin = zmin,
              zmax = zmax);
        
    if (mode == TRTT_PARALLEL_BEAM) {
        h_evaluator, p, "_trtt_single_spline_driven_PB2D_projector";
    }
        
    // else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM){
    //     h_evaluator, p, _trtt_single_spline_driven_FB2D_projector;
    // }

    else {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    return p;
}
F_single_2D = trtt_single_spline_driven_2D_projector;
// Define an easier alias

func _trtt_single_spline_driven_PB2D_projector(p, x, job)
/* DOCUMENT
*/           
{
    /* Copy scalars and make fast aliases for arrays. */
    nx = p.X.nx;
    ny = p.X.ny;
    s_scl = p.X.s_scl;
    s_deg = p.X.s_deg;
    nv = p.Yk.nv;
    v_scl = p.Yk.v_scl;
    theta = p.Yk.theta;
    local xos_coefs; eq_nocopy, xos_coefs, p.xos_coefs;
    local xsd_coefs; eq_nocopy, xsd_coefs, p.xsd_coefs;

    local z; eq_nocopy, z, p.z;        step_z = p.step_z;
    zmax = p.zmax;                     zmin = p.zmin;
    local w; eq_nocopy, w, p.w;        step_w= p.step_w;

    Np = p.X.size_footprint;
    local fft_D; eq_nocopy, fft_D, p.fft_D;

    /* Calculated parameters */
    xn = indgen(nx)-1.0;
    yn = indgen(ny)-1.0;
    v = trtt_pix_coord(p.Yk).v;
    costh = cos(theta);
    sinth = sin(theta);

    /* Allocate array for the result */
    if (!job) {
        out = array(double, nv);
    } else if (job == 1) {
        // "Object-dimensions-like" array
        out = array(double, nx, ny, 1); // FIXME: nt dimension equals
                                        // 1 but must be specified for
                                        // compatibility with time
                                        // interpolator
    } else {
        error, "unsupported JOB";
    }

    /* Fourier transform of the profile */
    // fft_R = (s_scl^2.0)*(sinc(s_scl*w))^(s_deg+1.0);
    /* Fourier transform of the profile */
    fft_R = (s_scl^2.0)*(sinc((s_scl*costh)*w)*sinc((s_scl*sinth)*w))^(s_deg+1.0);
    /* Multiply FFTs = convolution in Fourier space */
    fft_RD = complex(fft_R)*fft_D;
    /* Inverse FFT of the convolved Fourier profile */
    RD = step_w*double(roll(fft(complex(fft_RD), -1),(Np/2)-!(Np%2)));

    /* Project the object samples on the detector */
    xy_proj = _trtt_2D_projection_O_to_D_spline_driven (xn, yn, xos_coefs, xsd_coefs, TRTT_PARALLEL_BEAM);
    /* Position of the projected sample relative to the profile */
    pu = v(1) - xy_proj;
    
    /* FIXME: The profile is the same for each sample for a given projection
     * angle */
    
    for (j=1; j<=nv; ++j) {
        /* Evaluate which samples have an influence on the current pixel
         * detector */
        idp = where(pu>=zmin & pu<=zmax);
        if(is_array(idp)) {
            // idp = idp(sort(idp));
            /* Evaluate the corresponding index on the profile */
            idRD = long(floor(pu(idp)/step_z) + (Np/2) - !(Np%2) + 1);
            if (!job) {
                /* Direct operator : sum of the profile values at "hit"
                 * positions */
                out(j) = sum(x(idp,1)*RD(idRD));
            } else {
                /* Transpose operator */
                out(idp,1) += x(j)*RD(idRD);
            }
        }
        /* Increment the detector position */
        pu += v_scl;
    }
    
    /* Return data */
    return out;
}

func trtt_whole_spline_driven_2D_projector(X, Y, mode, status=)
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

        Ck = F_single_2D(X, Yk, mode, status=status);
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
F_whole_2D = trtt_whole_spline_driven_2D_projector;
