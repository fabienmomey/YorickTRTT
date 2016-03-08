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

func trtt_sparse_mult(x, MCk_coefs, MCk_row_dimlist, MCk_row_indices, MCk_col_dimlist, MCk_col_indices, job) {
    number = numberof(MCk_coefs);
    if (!job) {
        out = array(double, MCk_row_dimlist);
        trtt_mvmult, out, x, number, MCk_coefs, MCk_row_indices, MCk_col_indices;
    } else {
        out = array(double, MCk_col_dimlist);
        trtt_mvmult, out, x, number, MCk_coefs, MCk_col_indices, MCk_row_indices;
    }   
    return out;
}

func trtt_2D_mpy_eval(x, &gx, Cops)
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    
    /* Extraction des paramètres communs */
    C1 = h_get(Cops,Ck_list(1));
    X = C1.X;
    nx = X.nx; ny = X.ny;
    s_scl = X.s_scl;
    s_deg = X.s_deg;
    nv = C1.Yk.nv;
    v_scl = C1.Yk.v_scl;
    local footprint; eq_nocopy, footprint, X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, X.cum_footprint;
    size_footprint = X.size_footprint;
    step_footprint = X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, X.coord_footprint;

    /* Extraction de la fonction de projection */
    funcproj = name_of_symlink(C1.projector);

    /* Envoi des paramètres communs */
    mp_exec, "mp_handout, x, nx, ny, s_scl, s_deg, nv, v_scl, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, funcproj; mp_handin;";
    /* Récupération par tous les rangs de la fonction de projection */
    mp_exec, "__projector = symlink_to_name(funcproj); mp_handin;";

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
            mp_exec, "mp_handout, k, nproc; mp_handin;";
            /* Envoi des paramètres propres */
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(xos_coefs), vpack(xsd_coefs), vpack(weight), vpack(yk);} else if (mp_rank==k) {vunpack, mp_recv(0), xos_coefs; vunpack, mp_recv(0), xsd_coefs; vunpack, mp_recv(0), weight; vunpack, mp_recv(0), yk;} mp_handin;";
        }

        /* Paramètres propres du rang 0 */
        xos_coefs = Ck_0.xos_coefs;
        xsd_coefs = Ck_0.xsd_coefs;
        weight = Ck_0.Yk.Wk.a;
        yk = Ck_0.Yk.pixels;

        mp_exec, "if (mp_rank<nproc) {Cxk = array(double, nv); __projector, Cxk, x, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, 0; Cxmyk = Cxk-yk; WCxmyk = 0.5*sum(Cxmyk*weight*Cxmyk);} else {WCxmyk=0.0;} mp_handin;";
        mp_exec, "if (mp_rank<nproc) {WCxk = weight*Cxk; WCxkt = array(double, nx, ny); __projector, WCxkt, WCxk, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, 1;} else {WCxkt = 0.0;} mp_handin;";
        mp_exec, "residuals = mp_handin(WCxmyk); gradient = mp_handin(WCxkt(*)); mp_handin;";
        /* Residuals */
        fx += residuals;
        /* Gradient */
        gx += gradient;
    }
    
    return fx;
}

func trtt_2D_mpy_eval_sparse(x, &gx, Cops)
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    
    /* Extraction des paramètres communs */
    C1 = h_get(Cops,Ck_list(1));
    X = C1.X;
    nv = C1.Yk.nv;
    v_scl = C1.Yk.v_scl;

    /* Envoi des paramètres communs */
    mp_exec, "mp_handout, x; mp_handin;";

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
            weight = Ck.Yk.Wk.a;
            yk = Ck.Yk.pixels;
            MCk_coefs = Ck.MCk_coefs;
            MCk_row_dimlist = Ck.MCk_row_dimlist;
            MCk_row_indices = Ck.MCk_row_indices;
            MCk_col_dimlist = Ck.MCk_col_dimlist;
            MCk_col_indices = Ck.MCk_col_indices;

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc; mp_handin;";
            /* Envoi des paramètres propres */
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(MCk_coefs), vpack(MCk_row_dimlist), vpack(MCk_row_indices), vpack(MCk_col_dimlist), vpack(MCk_col_indices), vpack(weight), vpack(yk);} else if (mp_rank==k) {vunpack, mp_recv(0), MCk_coefs; vunpack, mp_recv(0), MCk_row_dimlist; vunpack, mp_recv(0), MCk_row_indices; vunpack, mp_recv(0), MCk_col_dimlist; vunpack, mp_recv(0), MCk_col_indices; vunpack, mp_recv(0), weight; vunpack, mp_recv(0), yk;} mp_handin;";
        }

        /* Paramètres propres du rang 0 */
        xos_coefs = Ck_0.xos_coefs;
        xsd_coefs = Ck_0.xsd_coefs;
        weight = Ck_0.Yk.Wk.a;
        yk = Ck_0.Yk.pixels;
        MCk_coefs = Ck_0.MCk_coefs;
        MCk_row_dimlist = Ck_0.MCk_row_dimlist;
        MCk_row_indices = Ck_0.MCk_row_indices;
        MCk_col_dimlist = Ck_0.MCk_col_dimlist;
        MCk_col_indices = Ck_0.MCk_col_indices;


        mp_exec, "if (mp_rank<nproc) {Cxk = trtt_sparse_mult(x, MCk_coefs, MCk_row_dimlist, MCk_row_indices, MCk_col_dimlist, MCk_col_indices, 0n); Cxmyk = Cxk-yk; WCxmyk = 0.5*sum(Cxmyk*weight*Cxmyk);} else {WCxmyk=0.0;} mp_handin;";
        mp_exec, "if (mp_rank<nproc) {WCxk = weight*Cxk; WCxkt = trtt_sparse_mult(WCxk, MCk_coefs, MCk_row_dimlist, MCk_row_indices, MCk_col_dimlist, MCk_col_indices, 1n);} else {WCxkt = 0.0;} mp_handin;";
        mp_exec, "residuals = mp_handin(WCxmyk); gradient = mp_handin(WCxkt(*)); mp_handin;";
        /* Residuals */
        fx += residuals;
        /* Gradient */
        gx += gradient;
    }
    
    return fx;
}

/* USER FUNCTIONS ============================================================ */

func trtt_2D_optim_cost_function_quadratic_mpy(x, &gx, cost_mat)
/* DOCUMENT trtt_optim_cost_function_quadratic, x, &gx, cost_mat
   
   This function calculates the new cost function FX, and its gradient GX, from
   an updated estimate X of the iterative minimization algorithm, in the context
   of 2D tomographic reconstruction using a set of independent projectors COPS
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

   The hash table COST_MAT contains the set of projectors COPS, which contains the
   data Y, the weight matrix W and eventually a regularization term (constraint
   function) REGUL.  If defined, the constraint function coming from REGUL is
   added to the cost function FX and its gradient to GX.
   
   SEE ALSO: trtt_optim_define_cost_mat.
*/
{
    dims = dimsof(x);
    /* Projector */
    Cops = cost_mat.Cops;

    WCyT = cost_mat.WCyT;
    
    /* Regularization */
    regul = cost_mat.regul;
    regulTV = cost_mat.regulTV;
    /* FIXME: covariance matrix not implemented yet */

    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        wx = cost_mat.wx;
        wy = cost_mat.wy;
        wz = cost_mat.wz;
        wt = cost_mat.wt;
    }

    gx = 0.0;
    if (!exist_sparse_coefs) {
        fx = trtt_2D_mpy_eval(x, gx, Cops);
    } else {
        fx = trtt_2D_mpy_eval_sparse(x, gx, Cops);
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
   
    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 2) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, 1);
            } else if (dims(1) == 3) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wt);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wt, 1);
            } else if (dims(1) == 4) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wz, wt);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, wt, 1);
            }
            
            
        } else {
            /* Adding to the global cost function */
            fx_plus = rgl_totvar(x, gx_plus, weight=weight, threshold=threshold, options=options);
        }

        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    }
    
    else if (!is_void(regul)) {

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 2) {
                z = spl_interp_apply(x, wx, wy);
                /* Updating the regularization term and clearing cash */
                rgl_update, regul, z;

                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);                               
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, 1);

                
            } else if (dims(1) == 3) {
                z = spl_interp_apply(x, wx, wy, wt);
                /* Updating the regularization term and clearing cash */
                rgl_update, regul, z;

                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);                               
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wt, 1);

                
            } else if (dims(1) == 4) {
                z = spl_interp_apply(x, wx, wy, wz, wt);
                /* Updating the regularization term and clearing cash */
                rgl_update, regul, z;

                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, wt, 1);
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
trtt_cost_quadratic_mpy = trtt_2D_optim_cost_function_quadratic_mpy; // Define an easier alias

func trtt_2D_optim_cost_function_quadratic_mpy_opky(cost_mat, x, &gx)
/* DOCUMENT trtt_optim_cost_function_quadratic, x, &gx, cost_mat
   
   This function calculates the new cost function FX, and its gradient GX, from
   an updated estimate X of the iterative minimization algorithm, in the context
   of 2D tomographic reconstruction using a set of independent projectors COPS
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

   The hash table COST_MAT contains the set of projectors COPS, which contains the
   data Y, the weight matrix W and eventually a regularization term (constraint
   function) REGUL.  If defined, the constraint function coming from REGUL is
   added to the cost function FX and its gradient to GX.
   
   SEE ALSO: trtt_optim_define_cost_mat.
*/
{
    dims = dimsof(x);
    /* Projector */
    Cops = cost_mat.Cops;

    WCyT = cost_mat.WCyT;
    
    /* Regularization */
    regul = cost_mat.regul;
    regulTV = cost_mat.regulTV;
    /* FIXME: covariance matrix not implemented yet */

    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        wx = cost_mat.wx;
        wy = cost_mat.wy;
        wz = cost_mat.wz;
        wt = cost_mat.wt;
    }

    gx = 0.0;
    if (!exist_sparse_coefs) {
        fx = trtt_2D_mpy_eval(x, gx, Cops);
    } else {
        fx = trtt_2D_mpy_eval_sparse(x, gx, Cops);
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
   
    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 2) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, 1);
            } else if (dims(1) == 3) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wt);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wt, 1);
            } else if (dims(1) == 4) {              
                /* Adding to the global cost function */
                z = spl_interp_apply(x, wx, wy, wz, wt);
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, wt, 1);
            }
            
            
        } else {
            /* Adding to the global cost function */
            fx_plus = rgl_totvar(x, gx_plus, weight=weight, threshold=threshold, options=options);
        }

        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    }
    
    else if (!is_void(regul)) {

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 2) {
                z = spl_interp_apply(x, wx, wy);
                /* Updating the regularization term and clearing cash */
                rgl_update, regul, z;

                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);                               
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, 1);

                
            } else if (dims(1) == 3) {
                z = spl_interp_apply(x, wx, wy, wt);
                /* Updating the regularization term and clearing cash */
                rgl_update, regul, z;

                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);                               
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wt, 1);

                
            } else if (dims(1) == 4) {
                z = spl_interp_apply(x, wx, wy, wz, wt);
                /* Updating the regularization term and clearing cash */
                rgl_update, regul, z;

                /* Adding to the global cost function */
                fx_plus = rgl_get_penalty(regul, z);
                
                /* Adding the gradient of the constraint to the gradient of the
                   cost function */
                gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, wt, 1);
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
trtt_cost_quadratic_mpy_opky = trtt_2D_optim_cost_function_quadratic_mpy_opky; // Define an easier alias

func trtt_2D_optim_define_cost_mat(Cops, dweights, use_sparse_coefs=, matrix=, regul=, regulTV=, viewer=, win_viewer=, win_viewer2=, cmin=, cmax=, roi=, xref=)
/* DOCUMENT CT_system_OptimPack_define_cost_mat, Cops, dweights, use_sparse_coefs=, matrix=, regul=, regulTV=, viewer=, win_viewer=, cmin=, cmax=, roi=, xref=
   
   This function defines the parameters for calculating the cost function from
   an updated estimate X of an iterative minimization algorithm for 2D
   tomographic reconstruction.

   The hash table COST_MAT contains the set of projectors COPS, which contains the
   data Y, the weight matrix W and eventually a regularization term (constraint
   function) REGUL.  REGUL is defined with the function RGL_NEW and
   RGL_CONFIG. If MATRIX=1, and if it is not done yet, the whole sparse
   projection matrix MCOPS is calculated, when not-on-the-fly projection
   calculation is asked. ROI and XREF are parameters used for metrics
   calculation. ROI defines a spatial region of interest, thus it is a 4 element
   array corresponding to the bounds (indexes) of the ROI :
   [id_min,idx_max,idy_min,idy_max].XREF is s reference image used metric RMSE
   is asked. VIEWER specifies the viewer function to apply at each iteration.
   
   SEE ALSO: trtt_2D_optim_simu_launcher, trtt_optim_cost_function_quadratic,
             rgl_new, rgl_config.
*/
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    C1 = h_get(Cops,Ck_list(1));
    X = C1.X;
    nx = X.nx;
    ny = X.ny;
    nv = C1.Yk.nv;

    cost_mat = h_new(
        Cops=Cops,
        regul=regul,
        regulTV=regulTV,
        viewer=viewer,
        win_viewer=win_viewer,
        win_viewer2=win_viewer2,
        cmin=cmin,
        cmax=cmax
        );
    
    if (matrix || Cops.matrix) {
        /* Get sparse matrix */
        MCops = _trtt_2D_sparse_matrix(Cops);
        /* Get sinogram data */
        y = [];
        weights = [];
        for (k=1; k<=ndata; ++k) {
            /* Get the current projector */
            Ck = h_get(Cops,Ck_list(k));
            /* Get and store the current projection data */
            Yk = Ck.Yk;
            grow, y, Yk.pixels;
            grow, weights, Yk.weight;
        }
        /* Set sparse matrix parameters in COST_MAT */
        h_set, cost_mat, MCops = MCops;
        h_set, cost_mat, y = y;
        if (!is_void(weights))
            h_set, cost_mat, W = linop_new(LINOP_DIAGONAL, dweights(*)*weights);
        else
            h_set, cost_mat, W = linop_new(LINOP_DIAGONAL, dweights(*)*array(1.0,numberof(y)));   
    } else {
        /* Calculate true data transposition */
        WCyT = 0.0;
        if ((!is_void(use_sparse_coefs)||use_sparse_coefs) & !Cops.exist_sparse_coefs) {
            _trtt_2D_sparse_coefs, Cops;
        }
        
        for (k=1; k<=ndata; ++k) {
            Ck = h_get(Cops, Ck_list(k));
            if (Cops.exist_sparse_coefs) {
                MCk = sparse_matrix(Ck.MCk_coefs, Ck.MCk_row_dimlist, Ck.MCk_row_indices, Ck.MCk_col_dimlist, Ck.MCk_col_indices);
                h_set, Ck, MCk=MCk;
            }
            Yk = Ck.Yk;
            yk = Yk.pixels;
            
            if (!Cops.roi_cov) {
                weight = Yk.weight;
                if (!is_void(weight)) {
                    Wk = linop_new(LINOP_DIAGONAL, dweights(,k)*weight);
                } else {
                    Wk = linop_new(LINOP_DIAGONAL, dweights(,k)*array(1.0,numberof(yk)));
                }
                
                h_set, Yk, Wk = Wk;
            } else {
                Wk=Yk.Wk;
            }
            
            WCyT += Ck(Wk(yk),1);
        }
        /* Set data transposition in COST_MAT */
        h_set, cost_mat, WCyT = WCyT;
    }
    
    /* Set interpolation coefficients for SPLINE DRIVEN only */
    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        cx=(indgen(0:nx*16)/(nx*16-1.0))*(nx-1.0);
        cy=(indgen(0:ny*16)/(ny*16-1.0))*(ny-1.0);
        
        wx = spl_interp_coefs_new_spline(cx, nx, s_deg);
        wy = spl_interp_coefs_new_spline(cy, ny, s_deg);
        h_set, cost_mat, wx = wx;
        h_set, cost_mat, wy = wy;
        
        // h_set, cost_mat, wx = X.wx;
        // h_set, cost_mat, wy = X.wy;
        h_set, cost_mat, wz = X.wz;
        h_set, cost_mat, wt = X.wt;
    }
    
    /* Check for "iterative" metrics calculation */
    extern TRTT_2D_RESIDUALS_ITER;
    extern TRTT_2D_REGULTERM_ITER;
    extern TRTT_2D_RMSE_ITER;

    if (TRTT_2D_RESIDUALS_ITER || TRTT_2D_REGULTERM_ITER || TRTT_2D_RMSE_ITER) {
        extern rec_results;
        /* Get "iterative" metrics */
        if (TRTT_2D_RESIDUALS_ITER) {
            h_set, rec_results, residuals_iter=[];
        }
        if (TRTT_2D_REGULTERM_ITER) {
            h_set, rec_results, regulterm_iter=[];
        }
        if (TRTT_2D_RMSE_ITER) {
            h_set, rec_results, rmse_iter=[];
            h_set, cost_mat, xref=xref;
            h_set, cost_mat, s_deg=X.s_deg;
            h_set, cost_mat, t_deg=X.t_deg;
        }
        
        /* Update COST_MAT with metrics parameters */
        h_set, cost_mat, roi=roi;
    }
    
    return cost_mat;
}

func trtt_2D_optim_simu_launcher(Cops, cost, x=, method=, xmin=, xmax=, bnd=, mem=, factr=, frtol=, fatol=, pgtol=, maxiter=, maxeval=, dweights=, get_costs=, metrics=, roi=, xref=, use_sparse_coefs=, matrix=, regul=, regulTV=, viewer=, win_viewer=, win_viewer2=, cmin=, cmax=, verbose=, store=)        
/* DOCUMENT trtt_optim_simu_launcher, Cops, cost, x=, method=, mem=, dweights=, metrics=, roi=, xref=, use_sparse_coefs=, matrix=, regul=, regulTV=, xmin=, xmax=, maxiter=, viewer=, win_viewer=, cmin=, cmax=, verbose=, store=

   Reconstruction algorithm launcher for TRTT code. The L-BFGS optimization
   algorithm is used (OP_MNB with no METHOD/MEM specification). The hashtab COPS
   include all the needed parameters to create the COST_MAT argument for the cost
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
                  trtt_2D_metrics_info to know how to define this parameter.           
   ROI= :         Some metrics can be calculated on a spatial region of interest
                  (ROI). Thus ROI is a 4 element array corresponding to the
                  bounds (indexes) of the ROI : [id_min,idx_max,idy_min,idy_max].
   XREF= :        A reference image for metric RMSE (mandatory).
   MATRIX= :      Ask for sparse matrix calculation (if not done yet).  
   REGUL= :       Specify the REGULARIZATION term.
   XMIN, XMAX= :  Lower/upper bounds for  X. Must be  conformable with X.
                  For instance with XMIN=0, the non-negative solution will be
                  returned.
   MAXITER= :     Maximum number of iterations (default: no limits).
   
   VIEWER= :      User defined subroutine to call every iteration to display the
                  solution X. The subroutine will be called as: viewer, x,
                  cost_mat; where X is the current solution and COST_MAT is the value
                  of keyword COST_MAT (which to see). If the viewer uses Yorick
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

   SEE ALSO: op_mnb, trtt_2D_metrics_info.
*/
{
    Ck_list = Cops.Ck_list;
    /* Resulting TOMOBJ */
    X = h_get(Cops,Ck_list(1)).X;

    if (is_void(viewer)) viewer=0n;

    if (is_void(x)) {
        local x; eq_nocopy, x, X.voxels;
    }
    dims = dimsof(x);
    
    // /* FIXME: METRICS EVALUATION */
    // local rec_results; rec_results = h_new();   
    // // printer = _trtt_2D_metrics_for_lbfgsb_solve;
    // printer = _trtt_2D_metrics;
    // if (!is_void(metrics)) {
    //     flag_metrics = 1n;
    //      /* Create the results hashtab */
    //     local TRTT_2D_RESIDUALS; TRTT_2D_RESIDUALS = metrics(1);
    //     local TRTT_2D_RESIDUALS_ITER; TRTT_2D_RESIDUALS_ITER = metrics(2);
    //     local TRTT_2D_REGULTERM; TRTT_2D_REGULTERM = metrics(3);
    //     local TRTT_2D_REGULTERM_ITER; TRTT_2D_REGULTERM_ITER = metrics(4);
    //     local TRTT_2D_RMSE; TRTT_2D_RMSE = metrics(5);
    //     local TRTT_2D_RMSE_ITER; TRTT_2D_RMSE_ITER = metrics(6);
    //     if ((TRTT_2D_RMSE_ITER || TRTT_2D_RMSE) & is_void(xref)) {
    //         TRTT_2D_RMSE = 0n;
    //         TRTT_2D_RMSE_ITER = 0n;
    //         trtt_error_display_msg, code=TRTT_ERR_MISSING_PARAMETER, msg="No reference image for calculating RMSE. Nothing done.";
    //     }
    // } else {
    //     flag_metrics = 0n;
    //     local TRTT_2D_RESIDUALS; TRTT_2D_RESIDUALS = 0n;
    //     local TRTT_2D_RESIDUALS_ITER; TRTT_2D_RESIDUALS_ITER = 0n;
    //     local TRTT_2D_REGULTERM; TRTT_2D_REGULTERM = 0n;
    //     local TRTT_2D_REGULTERM_ITER; TRTT_2D_REGULTERM_ITER = 0n;
    //     local TRTT_2D_RMSE; TRTT_2D_RMSE = 0n;
    //     local TRTT_2D_RMSE_ITER; TRTT_2D_RMSE_ITER = 0n;
    // }

    // if (get_costs) {
    //     extern costs_iter;
    //     extern costs_eval;
    //     if (!is_void(costs_iter))
    //         costs_iter=costs_iter(1:-1);
    //     if (!is_void(costs_eval))
    //         costs_eval=costs_eval(1:-1);
    // }
    
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
    }
    
    /* Define COST_MAT */
    if (is_void(dweights)) {
        ndata = numberof(Ck_list);
        nv = h_get(Cops,Ck_list(1)).Yk.nv;
        dweights = array(1.0, nv, ndata);
    }

    /* Check multi-regularizations conflict */
    //FIXME: either REGUL or REGULTV must be set at a time
    if (!is_void(regul) & !is_void(regulTV)) {
        trtt_error_display_msg, msg="Either REGUL or REGULTV must be set at a time.";
        return TRTT_FAILURE;
    }

    /* FIXME: Define materials for h_evaluator of the cost function */
    cost_mat = trtt_2D_optim_define_cost_mat(Cops, dweights, use_sparse_coefs=use_sparse_coefs, matrix=matrix, regul=regul, regulTV=regulTV, viewer=viewer, win_viewer=win_viewer, win_viewer2=win_viewer2, cmin=cmin, cmax=cmax, roi=roi, xref=xref);

    /* FIXME: METRICS EVALUATION */
    local rec_results; rec_results = h_new();   
    if (!is_void(metrics)) {
        if (cost==trtt_cost_quadratic_mpy) {
            printer = _trtt_2D_metrics;
        } else if (cost==trtt_cost_quadratic_mpy_opky) {
            printer = _trtt_2D_metrics_for_opky;
        }
        flag_metrics = 1n;
        /* Create the results hashtab */
        local TRTT_2D_RESIDUALS; TRTT_2D_RESIDUALS = metrics(1);
        local TRTT_2D_RESIDUALS_ITER; TRTT_2D_RESIDUALS_ITER = metrics(2);
        local TRTT_2D_REGULTERM; TRTT_2D_REGULTERM = metrics(3);
        local TRTT_2D_REGULTERM_ITER; TRTT_2D_REGULTERM_ITER = metrics(4);
        local TRTT_2D_RMSE; TRTT_2D_RMSE = metrics(5);
        local TRTT_2D_RMSE_ITER; TRTT_2D_RMSE_ITER = metrics(6);
        if ((TRTT_2D_RMSE_ITER || TRTT_2D_RMSE) & is_void(xref)) {
            TRTT_2D_RMSE = 0n;
            TRTT_2D_RMSE_ITER = 0n;
            trtt_error_display_msg, code=TRTT_ERR_MISSING_PARAMETER, msg="No reference image for calculating RMSE. Nothing done.";
        }
    } else if (is_void(metrics) & viewer) {
        flag_metrics = 0n;
        local TRTT_2D_RESIDUALS; TRTT_2D_RESIDUALS = 0n;
        local TRTT_2D_RESIDUALS_ITER; TRTT_2D_RESIDUALS_ITER = 0n;
        local TRTT_2D_REGULTERM; TRTT_2D_REGULTERM = 0n;
        local TRTT_2D_REGULTERM_ITER; TRTT_2D_REGULTERM_ITER = 0n;
        local TRTT_2D_RMSE; TRTT_2D_RMSE = 0n;
        local TRTT_2D_RMSE_ITER; TRTT_2D_RMSE_ITER = 0n;
        if (cost==trtt_cost_quadratic_mpy) {
            printer=std_viewer;
        } else if (cost==trtt_cost_quadratic_mpy_opky) {;
            printer=std_viewer_opky;
        }
    }

    if (get_costs) {
        extern costs_iter;
        extern costs_eval;
        if (!is_void(costs_iter))
            costs_iter=costs_iter(1:-1);
        if (!is_void(costs_eval))
            costs_eval=costs_eval(1:-1);
    }

    /**********************************************************************************************************/
    /* FIXME: Launch reconstruction */
    if (cost==trtt_cost_quadratic_mpy) {
        mem_copy, &x, op_mnb(cost, x, extra=cost_mat, xmin=xmin, xmax=xmax, method=method, mem=mem, maxiter=maxiter, maxeval=maxeval, frtol=frtol, fatol=fatol, viewer=printer, verb=verbose);
    } else if (cost==trtt_cost_quadratic_mpy_opky) { 
        /* FIXME: ADAPTATION POUR NOUVEL OPTIMPACK OPKY */
        minimizer=closure("trtt_cost_quadratic_mpy_opky",cost_mat);
        mem_copy, &x, opk_minimize(minimizer, x, mem=mem, lower=xmin, upper=xmax, maxiter=maxiter, maxeval=maxeval, printer=printer, verb=verbose);
    }
    /**********************************************************************************************************/

    /* Get "post-reconstruction" metrics */
    if (TRTT_2D_RESIDUALS) {
        residuals = _trtt_2D_residuals(x, Cops);
        h_set, rec_results, residuals=residuals;
    }
    
    /* Get back to samples for SPLINE DRIVEN */
    if (type == TRTT_SPLINE_DRIVEN) {
        if (dims(1) == 2)
            xrec=spl_interp_apply(x, wx, wy);
        else if (dims(1) == 3)
            xrec=spl_interp_apply(x, wx, wy, wt);
        else
            xrec=spl_interp_apply(x, wx, wy, wz, wt);
    } else {
        xrec=x;
    }
    
    if (TRTT_2D_REGULTERM) {
        regulterm = _trtt_2D_regulterm(xrec, regul=regul, regulTV=regulTV);
        h_set, rec_results, regulterm=regulterm;
    }
    if (TRTT_2D_RMSE) {
        rmse = _trtt_2D_rmse(xrec, xref, X.s_deg, X.t_deg, roi=roi);
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
    

    /* Set reconstructed image and reconstruction parameters in results
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
    if (get_costs) {
        h_set, rec_results, costs_iter=costs_iter;
        h_set, rec_results, costs_eval=costs_eval;
    }
    
    return rec_results;  
}

/* OBSOLETE FUNCTIONS ======================================================== */

// func trtt_2D_optim_cost_function_quadratic_matrix(x, &gx, cost_mat)
// /* DOCUMENT trtt_optim_cost_function_quadratic_matrix, x, &gx, cost_mat
   
//    This function calculates the new cost function FX, and its gradient GX, from
//    an updated estimate X of the iterative minimization algorithm, in the context
//    of 2D tomographic reconstruction using a sparse projection matrix MCOPS and
//    sinogram data Y. The cost function (or objective function), also called data
//    attachment term is here the WEIGHTED LEAST SQUARES CRITERION:

//    FX = || Cops.x-y ||² = 0.5 * transpose(MCops.x-y) * W * (MCops.x-y)
                                
//    where W is the weight matrix (inverse noise covariance matrix).

//    The gradient is : GX = || MCops.x-y ||² = transpose(MCops) * W * (MCops.x-y).

//    The hash table COST_MAT contains the sparse projector MCops, which contains the
//    data Y, the weight matrix W and eventually a regularization term (constraint
//    function) REGUL.  If defined, the constraint function coming from REGUL is
//    added to the cost function FX and its gradient to GX.
   
//    SEE ALSO: trtt_optim_define_cost_mat.
// */
// {
//     dims = dimsof(x);
//     /* Regularization */
//     regul = cost_mat.regul;
//     /* Matrix */
//     local MCops; eq_nocopy, MCops, cost_mat.MCops;
//     /* Sinogram */
//     y = cost_mat.y(*);
//     /* Weight operator (Inverse noise covariance matrix) */
//     W = cost_mat.W;
    
//     type = Cops.type;
    
//     if (type == TRTT_SPLINE_DRIVEN) {
//         wx = cost_mat.wx;
//         wy = cost_mat.wy;
//         wz = cost_mat.wz;
//         wt = cost_mat.wt;
//     }

//     Mx = MCops(x(*));
//     Mxmy = Mx-y;
//     fx = 0.5*sum(Mxmy*Mxmy);
//     gx = MCops(Mxmy,1);
//     gx = reform(gx,dims);

//     /* Adding a constraint function for regularization */
//     if (!is_void(regul)) {

//         if (type == TRTT_SPLINE_DRIVEN) {        
//             /* For regularization we have to work in samples space, and not in
//              * spline coefficients space */
//             if (dims(1) == 3) {
//                 z = spl_interp_apply(x, wx, wy, wt);
//                 /* Updating the regularization term and clearing cash */
//                 rgl_update, regul, z;

//                 /* Adding to the global cost function */
//                 fx_plus = rgl_get_penalty(regul, z);                               
                
//                 /* Adding the gradient of the constraint to the gradient of the
//                    cost function */
//                 gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wt, 1);

                
//             } else if (dims(1) == 4) {
//                 z = spl_interp_apply(x, wx, wy, wz, wt);
//                 /* Updating the regularization term  and clearing cash */
//                 rgl_update, regul, z;

//                 /* Adding to the global cost function */
//                 fx_plus = rgl_get_penalty(regul, z);                               
                
//                 /* Adding the gradient of the constraint to the gradient of the
//                    cost function */
//                 gx += spl_interp_apply(rgl_get_gradient(regul, z), wx, wy, wz, wt, 1);
//             }

//         } else {
//             // Updating the regularization term
//             rgl_update, regul, x;
//             // Adding to the global cost function
//             fx_plus = rgl_get_penalty(regul, x);                               
            
//             // Adding the gradient of the constraint to the gradient of the
//             // cost function
//             gx += rgl_get_gradient(regul, x);
//         }
        
//         /* Adding to the global cost function */
//         fx += fx_plus;
//     }

//     /* Return the global cost function */
//     return fx;
// }
// trtt_cost_quadratic_matrix = trtt_2D_optim_cost_function_quadratic_matrix; // Define an easier alias
