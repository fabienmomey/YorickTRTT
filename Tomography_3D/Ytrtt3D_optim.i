/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt2D_optim.i --
 *
 * TRTT 2D optimization and reconstruction package.
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

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

func trtt_3D_mpy_eval(x, &gx, Cops)
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    
    /* Extraction des paramètres communs */
    C1 = h_get(Cops,Ck_list(1));
    X = C1.X;
    nx = X.nx; ny = X.ny; nz = X.nz;
    s_scl = X.s_scl;
    s_deg = X.s_deg;
    nv = C1.Yk.nv;
    nu = C1.Yk.nu;
    v_scl = C1.Yk.v_scl;
    u_scl = C1.Yk.u_scl;
    local footprint; eq_nocopy, footprint, X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, X.cum_footprint;
    size_footprint = X.size_footprint;
    step_footprint = X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, X.coord_footprint;

    /* Extraction de la fonction de projection */
    funcproj = name_of_symlink(C1.projector);

    /* Envoi des paramètres communs */
    mp_exec, "mp_handout, x, nx, ny, nz, s_scl, s_deg, nu, nv, u_scl, v_scl, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, funcproj;";
    /* Récupération par tous les rangs de la fonction de projection */
    mp_exec, "__projector = symlink_to_name(funcproj);";

    /* Calcul du nombre d'occurrences */
    Nmult = ndata/mp_size;
    Nrest = ndata%mp_size;

    /* Initialisation des résidus */
    fx = 0.0;

    /* Calcul des occurrences */
    for (i=0; i<=Nmult; ++i) {

        if (i<Nmult) nproc = mp_size;
        else if (Nrest > 0) nproc = Nrest;
        else break;
        
        k_0 = i*mp_size + 1;
        Ck_0 = h_get(Cops,Ck_list(k_0));

        for (k=1; k<nproc; ++k) {
            Ck = h_get(Cops,Ck_list(k_0+k));
            xos_coefs = Ck.xos_coefs;
            xsd_coefs = Ck.xsd_coefs;
            weight = Ck.Yk.Wk.a;
            yk = Ck.Yk.pixels;

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc";
            /* Envoi des paramètres propres */
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(xos_coefs), vpack(xsd_coefs), vpack(weight), vpack(yk);} else if (mp_rank==k) {vunpack, mp_recv(0), xos_coefs; vunpack, mp_recv(0), xsd_coefs; vunpack, mp_recv(0), weight; vunpack, mp_recv(0), yk;} ";
        }

        /* Paramètres propres du rang 0 */
        xos_coefs = Ck_0.xos_coefs;
        xsd_coefs = Ck_0.xsd_coefs;
        weight = Ck_0.Yk.Wk.a;
        yk = Ck_0.Yk.pixels;

        mp_exec, "if (mp_rank<nproc) {Cxk = array(double, nu, nv); __projector, Cxk, x, nx, ny, nz, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nu, nv, u_scl, v_scl, xos_coefs, xsd_coefs, 0; Cxmyk = Cxk-yk; WCxmyk = 0.5*sum(Cxmyk*weight*Cxmyk);} else {WCxmyk=0.0;}";
        
        mp_exec, "if (mp_rank<nproc) {WCxk = weight*Cxk; WCxkt = array(double, nx, ny, nz); __projector, WCxkt, WCxk, nx, ny, nz, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nu, nv, u_scl, v_scl, xos_coefs, xsd_coefs, 1;} else {WCxkt = 0.0;}";
        mp_exec, "residuals = mp_handin(WCxmyk); gradient = mp_handin(WCxkt(*));";
        /* Residuals */
        fx += residuals;
        /* Gradient */
        gx += gradient;
    }
    
    return fx;
}

/* USER FUNCTIONS ============================================================ */

func trtt_3D_optim_cost_function_quadratic(x, &gx, extra)
/* DOCUMENT trtt_optim_cost_function_quadratic, x, &gx, extra
   
   This function calculates the new cost function FX, and its gradient GX, from
   an updated estimate X of the iterative minimization algorithm, in the context
   of 3D tomographic reconstruction using a set of independent projectors COPS
   and sinogram data Y. The cost function (or objective function), also called
   data attachment term is here the WEIGHTED LEAST SQUARES CRITERION:

   FX = || Cops.x-y ||² = 0.5 * transpose(Cops.x-y) * W * (Cops.x-y)
                                ___
                                \
                        = 0.5 * /__ || Ck.xk-yk ||²

                                ___
                                \
                        = 0.5 * /__ transpose(Ck.xk-yk) * Wk * (Ck.xk-yk)
                        
                                
   where W is the weight matrix (inverse noise covariance matrix).

   The gradient is : GX = || Cops.x-y ||² = transpose(Cops) * W * (Cops.x-y).

                                ___
                                \
                        =       /__ transpose(Ck) * Wk * (Ck.xk-yk)

   The hash table EXTRA contains the set of projectors COPS, which contains the
   data Y, the weight matrix W and eventually a regularization term (constraint
   function) REGUL.  If defined, the constraint function coming from REGUL is
   added to the cost function FX and its gradient to GX.
   
   SEE ALSO: trtt_optim_define_extra.
*/
{
    /* Projector */
    Cops = extra.Cops;
    /* True data transposition */
    WCyT = extra.WCyT;
        
    /* Regularization */
    regul = extra.regul;
    regulTV = extra.regulTV;
    /* FIXME: covariance matrix not implemented yet */

    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        wx = extra.wx;
        wy = extra.wy;
        wz = extra.wz;
        wt = extra.wt;
    }

    Ck_list = h_get(Cops,"Ck_list");
    ndata = numberof(Ck_list);
    fx = 0.0;
    gx = 0.0;
    for (k=1; k<=ndata; ++k) {
        /* Get the projector */
        Ck = h_get(Cops, Ck_list(k));
        Yk = Ck.Yk;
        /* Get corresponding data */
        yk = Yk.pixels;
        /* Get noise covariance matrix */
        Wk = Yk.Wk;
        
        Cxk = Ck(x);
        Cxmyk = Cxk-yk;
        WCxmyk = Wk(Cxmyk);
        fx += 0.5*sum(Cxmyk*WCxmyk);
        gx += Ck(Wk(Cxk),1);
    }
    gx -= WCyT;

    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        lweight = regulTV.lweight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wz);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options, lweight=lweight);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, 1);
            } else if (dims(1) == 4) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wz, wt);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options, lweight=lweight);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, wt, 1);
            }

        } else {
            /* Adding to the global cost function */
            fx_plus = rgl_totvar(x, gx_plus, weight=weight, threshold=threshold, options=options, lweight=lweight);
        }

        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    } else if (!is_void(regul)) {
        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 4) {
                z = spl_interp_apply(x, wx, wy, wz, wt);
                /* Updating the regularization term  and clearing cash */
                rgl_update, regul, z;
                
                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, wt, 1);
            } else if (dims(1) == 3) {
                z = spl_interp_apply(x, wx, wy, wz);
                /* Updating the regularization term  and clearing cash */
                rgl_update, regul, z;
                
                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, 1);
            }
        } else {
            // Updating the regularization term
            rgl_update, regul, x;
            // Adding to the global cost function
            fx_plus = rgl_get_penalty(regul, x);
            
            // Adding the gradient of the constraint to the gradient of the
            // cost function
            gx += rgl_get_gradient(regul, x);
        }

        /* Adding to the global cost function */
        fx += fx_plus;
    }

    /* Return the global cost function */
    return fx;
}
trtt3D_cost_quadratic = trtt_3D_optim_cost_function_quadratic; // Define an easier alias

func trtt_3D_optim_cost_function_quadratic_mpy(x, &gx, extra)
/* DOCUMENT trtt_optim_cost_function_quadratic_mpy, x, &gx, extra
   
   This function calculates the new cost function FX, and its gradient GX, from
   an updated estimate X of the iterative minimization algorithm, in the context
   of 3D tomographic reconstruction using a set of independent projectors COPS
   and sinogram data Y. The cost function (or objective function), also called
   data attachment term is here the WEIGHTED LEAST SQUARES CRITERION:

   FX = || Cops.x-y ||² = 0.5 * transpose(Cops.x-y) * W * (Cops.x-y)
                                ___
                                \
                        = 0.5 * /__ || Ck.xk-yk ||²

                                ___
                                \
                        = 0.5 * /__ transpose(Ck.xk-yk) * Wk * (Ck.xk-yk)
                        
                                
   where W is the weight matrix (inverse noise covariance matrix).

   The gradient is : GX = || Cops.x-y ||² = transpose(Cops) * W * (Cops.x-y).

                                ___
                                \
                        =       /__ transpose(Ck) * Wk * (Ck.xk-yk)

   The hash table EXTRA contains the set of projectors COPS, which contains the
   data Y, the weight matrix W and eventually a regularization term (constraint
   function) REGUL.  If defined, the constraint function coming from REGUL is
   added to the cost function FX and its gradient to GX.
   
   SEE ALSO: trtt_optim_define_extra.
*/
{
    dims = dimsof(x);
    /* Projector */
    Cops = extra.Cops;
    /* True data transposition */
    WCyT = extra.WCyT;
        
    /* Regularization */
    regul = extra.regul;
    regulTV = extra.regulTV;
    /* FIXME: covariance matrix not implemented yet */

    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        wx = extra.wx;
        wy = extra.wy;
        wz = extra.wz;
        wt = extra.wt;
    }

    Ck_list = h_get(Cops,"Ck_list");
    ndata = numberof(Ck_list);
    
    gx = 0.0;
    fx = trtt_3D_mpy_eval(x, gx, Cops);
    gx -= WCyT(*);
    gx = reform(gx,dims);

    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        lweight = regulTV.lweight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wz);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options, lweight=lweight);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, 1);
            } else if (dims(1) == 4) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wz, wt);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options, lweight=lweight);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, wt, 1);
            }

        } else {
            /* Adding to the global cost function */
            fx_plus = rgl_totvar(x, gx_plus, weight=weight, threshold=threshold, options=options, lweight=lweight);
        }

        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    } else if (!is_void(regul)) {
        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 4) {
                z = spl_interp_apply(x, wx, wy, wz, wt);
                /* Updating the regularization term  and clearing cash */
                rgl_update, regul, z;
                
                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, wt, 1);
            } else if (dims(1) == 3) {
                z = spl_interp_apply(x, wx, wy, wz);
                /* Updating the regularization term  and clearing cash */
                rgl_update, regul, z;
                
                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, 1);
            }
        } else {
            // Updating the regularization term
            rgl_update, regul, x;
            // Adding to the global cost function
            fx_plus = rgl_get_penalty(regul, x);
            
            // Adding the gradient of the constraint to the gradient of the
            // cost function
            gx += rgl_get_gradient(regul, x);
        }

        /* Adding to the global cost function */
        fx += fx_plus;
    }

    /* Return the global cost function */
    return fx;
}
trtt3D_cost_quadratic_mpy = trtt_3D_optim_cost_function_quadratic_mpy; // Define an easier alias

func trtt_3D_optim_define_extra(Cops, dweights, regul=, regulTV=, viewer=, win_viewer=, cmin=, cmax=, roi=, xref=)
/* DOCUMENT CT_system_OptimPack_define_extra, Cops, dweights, rec_results, roi=, regul=, regulTV=, viewer=, win_viewer=, cmin=, cmax=, roi=, xref=
   
   This function defines the parameters for calculating the cost function from
   an updated estimate X of an iterative minimization algorithm for 3D
   tomographic reconstruction.

   The hash table EXTRA contains the set of projectors COPS, which contains the
   data Y, the weight matrix W and eventually a regularization term (constraint
   function) REGUL.  REGUL is defined with the function RGL_NEW and
   RGL_CONFIG. If MATRIX=1, and if it is not done yet, the whole sparse
   projection matrix MCOPS is calculated, when not-on-the-fly projection
   calculation is asked. ROI and XREF are parameters used for metrics
   calculation. ROI defines a spatial region of interest, thus it is a 4 element
   array corresponding to the bounds (indexes) of the ROI :
   [id_min,idx_max,idy_min,idy_max].XREF is s reference image used metric RMSE
   is asked. VIEWER specifies the viewer function to apply at each iteration.
   
   SEE ALSO: trtt_3D_optim_simu_launcher, trtt_optim_cost_function_quadratic,
             rgl_new, rgl_config.
*/
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    X = h_get(Cops,Ck_list(1)).X;
    nx = X.nx;
    ny = X.ny;
    nz = X.nz;
    nu = h_get(Cops,Ck_list(1)).Yk.nu;
    nv = h_get(Cops,Ck_list(1)).Yk.nv;

    extra = h_new(
        Cops=Cops,
        regul=regul,
        regulTV=regulTV,
        nz=nz,
        viewer=viewer,
        win_viewer=win_viewer,
        cmin=cmin,
        cmax=cmax
        );
    
    /* Calculate true data transposition */
    WCyT = 0.0;
    for (k=1; k<=ndata; ++k) {
        Ck = h_get(Cops, Ck_list(k));
        Yk = Ck.Yk;
        nu = Yk.nu;
        nv = Yk.nv;
        yk = Yk.pixels;
        weight = Yk.weight;
        if (!is_void(weight))
            Wk = linop_new(LINOP_DIAGONAL, dweights(..,k)*weight);
        else
            Wk = linop_new(LINOP_DIAGONAL, dweights(..,k)*array(1.0,nu,nv));
        
        WCyT += Ck(Wk(yk),1);
        
        h_set, Yk, Wk = Wk;
    }
    /* Set data transposition in EXTRA */
    h_set, extra, WCyT = WCyT;
    
    // FIXME: Data covariance matrix framework not yet implemented

    /* Set interpolation coefficients for SPLINE DRIVEN only */
    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        h_set, extra, wx = X.wx;
        h_set, extra, wy = X.wy;
        h_set, extra, wz = X.wz;
        h_set, extra, wt = X.wt;
    }
    
    /* Check for "iterative" metrics calculation */
    extern TRTT_3D_RESIDUALS_ITER;
    extern TRTT_3D_REGULTERM_ITER;
    extern TRTT_3D_RMSE_ITER;

    if (TRTT_3D_RESIDUALS_ITER || TRTT_3D_REGULTERM_ITER || TRTT_3D_RMSE_ITER) {
        extern rec_results;
        /* Get "iterative" metrics */
        if (TRTT_3D_RESIDUALS_ITER) {
            h_set, rec_results, residuals_iter=[];
        }
        if (TRTT_3D_REGULTERM_ITER) {
            h_set, rec_results, regulterm_iter=[];
        }
        if (TRTT_3D_RMSE_ITER) {
            h_set, rec_results, rmse_iter=[];
            h_set, extra, xref=xref;
            h_set, extra, s_deg=X.s_deg;
            h_set, extra, t_deg=X.t_deg;
        }
        
        /* Update EXTRA with metrics parameters */
        h_set, extra, roi=roi;
    }
    
    return extra;
}

func trtt_3D_optim_simu_launcher(Cops, cost, x=, method=, mem=, dweights=, metrics=, roi=, xref=, weight=, regul=, regulTV=, xmin=, xmax=, maxiter=, viewer=, win_viewer=, cmin=, cmax=, verbose=, store=)        
/* DOCUMENT trtt_3D_optim_simu_launcher, Cops, cost, x=, method=, mem=, dweights=, metrics=, roi=, xref=, weight=, regul=, regulTV=, xmin=, xmax=, maxiter=, viewer=, win_viewer=, cmin=, cmax=, verbose=, store=   

   Reconstruction algorithm launcher for TRTT code. The L-BFGS optimization
   algorithm is used (OP_MNB with no METHOD/MEM specification). The hashtab COPS
   include all the needed parameters to create the EXTRA argument for the cost
   function COST: the projectors Ck, the data objects Yk and the resulting
   tomobj X. Other parameters can be specified to be passed to OP_MNB
   optimization algorithm.
   METHOD= :      Scalar integer which defines the optimization method to use.
                  Conjugate gradient algorithm is used if one of the bits
                  OP_FLAG_POLAK_RIBIERE, OP_FLAG_FLETCHER_REEVES, or
                  OP_FLAG_HESTENES_STIEFEL is set; otherwise, a limited memory variable
                  metric algorithm (VMLM-B) is used.  If METHOD is not specified and if
                  MEM=0, a conjugate gradient search is attempted with flags:
                  (OP_FLAG_UPDATE_WITH_GP |
                   OP_FLAG_SHANNO_PHUA    |
                   OP_FLAG_MORE_THUENTE   |
                   OP_FLAG_POLAK_RIBIERE  |
                   OP_FLAG_POWELL_RESTART)           
                  otherwise VMLM-B is used with flags:
                  (OP_FLAG_UPDATE_WITH_GP |
                   OP_FLAG_SHANNO_PHUA    |
                   OP_FLAG_MORE_THUENTE).                 
                  See documentation of op_get_flags to figure out the allowed
                  bit flags and their meaning.
   MEM= :         Number of previous directions used in variable metric limited memory
                  method (default min(7, numberof(X))).  
   METRICS= :     Table specifying the evaluation metrics to calculate (it is
                  possible not to calculate anyone) : see the help of
                  trtt_3D_metrics_info to know how to define this parameter.           
   ROI= :         Some metrics can be calculated on a spatial region of interest
                  (ROI). Thus ROI is a 4 element array corresponding to the
                  bounds (indexes) of the ROI : [id_min,idx_max,idy_min,idy_max].
   XREF= :        A reference image for metric RMSE (mandatory). 
   WEIGHT= :      (Not available) Specify the WEIGHT for EXTRA argument of the cost
                  function.
   REGUL= :       Specify the REGULARIZATION term.
   XMIN, XMAX= :  Lower/upper bounds for  X. Must be  conformable with X.
                  For instance with XMIN=0, the non-negative solution will be
                  returned.
   MAXITER= :     Maximum number of iterations (default: no limits).
   
   VIEWER= :      User defined subroutine to call every iteration to display the
                  solution X. The subroutine will be called as: viewer, x,
                  extra; where X is the current solution and EXTRA is the value
                  of keyword EXTRA (which to see). If the viewer uses Yorick
                  graphics window(s) it may call "pause, 1;" before returning to
                  make sure that graphics get correctly updated.
             
   WIN_VIEWER= :  Number of the display window. Default is 0.
   CMIN= :        minimum value for display. Default is min(estimate).
   CMAX= :        maximum value for display. Default is max(estimate).
   
   VERBOSE= :     Verbose mode?  If non-nil and non-zero, print out information
                  every VERB iterations and for the final one.
   STORE= :       Store reconstruction image in X.

   Depending on the cost function COST, the used projectors can be functions
   which perform the voxels projections on-the-fly, or sparse matrix operators
   which stores the pre-calculated coefficients of voxels projection.

   SEE ALSO: op_mnb, trtt_3D_metrics_info.
*/
{
    Ck_list = Cops.Ck_list;
    /* Resulting TOMOBJ */
    X = h_get(Cops,Ck_list(1)).X;

    if (is_void(x)) {
        local x; eq_nocopy, x, X.voxels;
    }
    dims = dimsof(x);

    /* FIXME: METRICS EVALUATION */
    local rec_results; rec_results = h_new();
    markup = _trtt_3D_metrics;
    if (!is_void(metrics)) {
        flag_metrics = 1n;
         /* Create the results hashtab */
        local TRTT_3D_RESIDUALS; TRTT_3D_RESIDUALS = metrics(1);
        local TRTT_3D_RESIDUALS_ITER; TRTT_3D_RESIDUALS_ITER = metrics(2);
        local TRTT_3D_REGULTERM; TRTT_3D_REGULTERM = metrics(3);
        local TRTT_3D_REGULTERM_ITER; TRTT_3D_REGULTERM_ITER = metrics(4);
        local TRTT_3D_RMSE; TRTT_3D_RMSE = metrics(5);
        local TRTT_3D_RMSE_ITER; TRTT_3D_RMSE_ITER = metrics(6);
        if ((TRTT_3D_RMSE_ITER || TRTT_3D_RMSE) & is_void(xref)) {
            TRTT_3D_RMSE = 0n;
            TRTT_3D_RMSE_ITER = 0n;
            trtt_error_display_msg, code=TRTT_ERR_MISSING_PARAMETER, msg="No reference image for calculating RMSE. Nothing done.";
        }
    } else {
        flag_metrics = 0n;
        local TRTT_3D_RESIDUALS; TRTT_3D_RESIDUALS = 0n;
        local TRTT_3D_RESIDUALS_ITER; TRTT_3D_RESIDUALS_ITER = 0n;
        local TRTT_3D_REGULTERM; TRTT_3D_REGULTERM = 0n;
        local TRTT_3D_REGULTERM_ITER; TRTT_3D_REGULTERM_ITER = 0n;
        local TRTT_3D_RMSE; TRTT_3D_RMSE = 0n;
        local TRTT_3D_RMSE_ITER; TRTT_3D_RMSE_ITER = 0n;
    }
    
    /* Get type of projector */
    type = Cops.type;
    /* Conversion to spline coefficients if necessary */
    if (type == TRTT_SPLINE_DRIVEN) {
        s_deg = X.s_deg;
        t_deg = X.t_deg;
        wx = X.wx;
        wy = X.wy;
        wz = X.wz;
        wt = X.wt;

        // if (dims(1) == 3)
        //     mem_copy, &x, spl_spline_samples_to_coefficients(x, s_deg, s_deg, s_deg);
        // else if (dims(1) == 4)
        //     mem_copy, &x, spl_spline_samples_to_coefficients(x, s_deg, s_deg, s_deg, t_deg);
    }
    
    /* Define EXTRA */
    if (is_void(dweights)) {
        ndata = numberof(Ck_list);
        nv = h_get(Cops,Ck_list(1)).Yk.nv;
        nu = h_get(Cops,Ck_list(1)).Yk.nu;
        dweights = array(1.0, nu, nv, ndata);
    }

    /* Check multi-regularizations conflict */
    //FIXME: either REGUL or REGULTV must be set at a time
    if (!is_void(regul) & !is_void(regulTV)) {
        trtt_error_display_msg, msg="Either REGUL or REGULTV must be set at a time.";
        return TRTT_FAILURE;
    }
    
    extra = trtt_3D_optim_define_extra(Cops, dweights, regul=regul, regulTV=regulTV, viewer=viewer, win_viewer=win_viewer, cmin=cmin, cmax=cmax, roi=roi, xref=xref);

    /* FIXME: Launch reconstruction */
    mem_copy, &x, op_mnb(cost, x, extra=extra, xmin=xmin, xmax=xmax, method=method, mem=mem, maxiter=maxiter, viewer=markup, verb=verbose);

    /* Get "post-reconstruction" metrics */
    if (TRTT_3D_RESIDUALS) {
        residuals = _trtt_3D_residuals(x, Cops);
        h_set, rec_results, residuals=residuals;
    }

    /* Get back to samples for SPLINE DRIVEN */
    if (type == TRTT_SPLINE_DRIVEN) {
        if (dims(1) == 3)
            xrec=spl_interp_apply(x, wx, wy, wz);
        else if (dims(1) == 4)
            xrec=spl_interp_apply(x, wx, wy, wz, wt);
    } else {
        xrec=x;
    }
    
    if (TRTT_3D_REGULTERM) {
        regulterm = _trtt_3D_regulterm(xrec, regul=regul, regulTV=regulTV);
        h_set, rec_results, regulterm=regulterm;
    }
    if (TRTT_3D_RMSE) {
        rmse = _trtt_3D_rmse(xrec, xref, roi=roi);
        h_set, rec_results, rmse=rmse;
    }
    
    /* Set reconstructed image in resulting TOMOBJ X */
    // if (!is_void(store)) trtt_tomobj_set_voxels, X, xr;

    /* Calculate mean noise standard deviation (for comparison with data
     * residuals) */
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    mean_noise_std = 0.0;
    for (k=1; k<=ndata; ++k) {
        Ck = h_get(Cops, Ck_list(k));
        mean_noise_std += Ck.Yk.mean_var_noise;
    }
    mean_noise_std /= double(ndata);
    mean_noise_std = sqrt(mean_noise_std);

    /* Set reconstructed image and recosntruction parameters in results
     * hashtab */
    h_set, rec_results, mean_noise_std = mean_noise_std;
    h_set, rec_results, regul=regul;
    h_set, rec_results, regulTV=regulTV;
    h_set, rec_results, maxiter=maxiter;
    h_set, rec_results, xmin=xmin;
    h_set, rec_results, xmax=xmax;
    h_set, rec_results, roi=roi;
    h_set, rec_results, xref=xref;
    h_set_copy, rec_results, "x", xrec;
    
    return rec_results;  
}

func trtt_3D_optim_standard_viewer(x, extra)
/* DOCUMENT trtt_optim_standard_viewer, x, extra
   
   Viewer for the call of the optimization algorithm for the reconstruction of a
   phantom. Plot the current estimate X.

*/
{
    nz = extra.nz;
    kz = long(floor(nz/2.0));
    nwin = extra.win_viewer;
    cmin = extra.cmin;
    cmax = extra.cmax;

    if (is_void(cmin)) cmin = min(x(*));
    if (is_void(cmax)) cmax = max(x(*));

    if (is_void(nwin)) nwin=0;
        
    mywindow, nwin; fma; img_plot, x(..,kz), cbar=1, cmin=cmin, cmax=cmax, vert=1, format="%.3f";

    pause, 1;
}
std_viewer = trtt_3D_optim_standard_viewer; // Define an easier alias

func trtt_3D_optim_low_contrast_viewer(x, extra)
/* DOCUMENT trtt_optim_low_constrast_viewer, x, extra
   
   Viewer for the call of the optimization algorithm for the reconstruction of a
   low contrast phantom (the image X will be histogram equalized), for example
   the Shepp Logan phantom. Plot the current estimate X.
*/
{
    nz = extra.nz;
    kz = long(floor(nz/2.0));
    nwin = extra.win_viewer;
    cmin = extra.cmin;
    cmax = extra.cmax;

    if (is_void(cmin)) cmin = min(x);
    if (is_void(cmax)) cmax = max(x);

    if (is_void(nwin)) nwin=0;
        
    trtt_3D_egalize_and_display_low_contrast, x(..,kz), cmin=cmin, cmax=cmax, nwin=nwin;

    pause, 1;
}
low_ctr_viewer = trtt_3D_optim_low_contrast_viewer; // Define an easier alias

func trtt_3D_egalize_and_display_low_contrast(z, first=, scale=, top=, cmin=, cmax=, nticks=, nwin=)
/* DOCUMENT trtt_egalize_and_display_low_contrast, z, top=, cmin=, cmax=, nticks=, nwin=
   
   This function is based on the HISTEQ_SCALE function which applies an
   egalization of the histogram of Z. Then the resulting image is plotted with
   the appropriate colorbar involving the "egalized" levels.  NTICKS= : Number
   of plotted labels on the colorbar.  NWIN= : Number of the display
   window. Default is 0.
   
   SEE ALSO: histeq_scale, img_plot.
*/
{
    if (is_void(nticks)) nticks = 11;
    if (is_void(nwin)) nwin=0;
    if (is_void(cmin)) cmin = min(z(*));
    if (is_void(cmax)) cmax= max(z(*));

        
    if (is_void(top)) top= bytscl([0.,1.])(2);  /* palette size - 1 */
    top= long(top);
    if (top<0 | top>255) error, "top value out of range 0-255";
    y= z(*);
    if (!is_void(cmin)) y= y(where(y>=cmin));
    if (!is_void(cmax)) y= y(where(y<=cmax));
    y= y(sort(y));
    x= span(0.,1., numberof(y));
    xp= span(0.,1., top+2);
    bins= interp(y, x, xp);
    list= where(bins(dif)<=0.0);
    if (numberof(list)) {
        /* some value (or values) of z are repeated many times --
           try to handle this by adding a small slope to the sorted y */
        dy= y(0)-y(1);
        if (!dy) dy= 1.0;
        for (eps=1.e-10 ; eps<1000.1 ; eps*=10.) {
            bins= interp(y+eps*dy*x, x, xp);
            list= where(bins(dif)<=0.0);
            if (!numberof(list)) break;
        }
        if (eps>1000.) error, "impossible error??";
    }

    // Apply the contrast egalization
    zs = char(max(min(digitize(z,bins)-2,top),0));

    // Define the labels
    levstep = int((top+2)/nticks)+1;
    labels = bins(1:0:levstep);
       
    // Plot the egalized image with the appropriate colorbar
    mywindow, nwin; fma; img_plot, zs, first=first, scale=scale, cbar=1, cmin=bins(1), cmax=bins(0), vert=1, nticks=nticks, format="%.3f", labels=labels;
}
trtt_pli_low_ctr = trtt_3D_egalize_and_display_low_contrast; // define an easier alias
