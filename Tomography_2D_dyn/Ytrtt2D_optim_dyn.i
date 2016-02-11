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

func trtt_2D_mpy_dyn_eval(x, &gx, Cops)
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    
    /* Extraction des paramètres communs */
    C1 = h_get(Cops,Ck_list(1));
    A1 = C1.Ak;
    B1 = C1.Bk;
    X = A1.X;
    nx = X.nx; ny = X.ny; nt = X.nt;
    s_scl = X.s_scl;
    s_deg = X.s_deg;
    nv = A1.Yk.nv;
    v_scl = A1.Yk.v_scl;
    local footprint; eq_nocopy, footprint, X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, X.cum_footprint;
    size_footprint = X.size_footprint;
    step_footprint = X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, X.coord_footprint;

    /* Extraction de la fonction de projection */
    funcproj = name_of_symlink(A1.projector);
    /* Extraction de la fonction d'interpolation temporelle */
    funcinterp = name_of_symlink(B1.interpolator);

    /* Envoi des paramètres communs */
    mp_exec, "mp_handout, nx, ny, nt, s_scl, s_deg, nv, v_scl, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, funcproj, funcinterp; mp_handin;";
    /* Récupération par tous les rangs de la fonction de projection */
    mp_exec, "__projector = symlink_to_name(funcproj); mp_handin;";
    mp_exec, "__interpolator = symlink_to_name(funcinterp); mp_handin;";

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
        Ak_0 = Ck_0.Ak;
        Bk_0 = Ck_0.Bk;
        p_0 = Bk_0.wt.p;
        j_0 = Bk_0.wt.j + 1;
        w_0 = Bk_0.wt.w;

        for (k=1; k<nproc; ++k) {
            Ck = h_get(Cops,Ck_list(k_0+k));
            /* Projector */
            Ak = Ck.Ak;
            xos_coefs = Ak.xos_coefs;
            xsd_coefs = Ak.xsd_coefs;
            weight = Ak.Yk.Wk.a;
            yk = Ak.Yk.pixels;
            /* Interpolator */
            Bk = Ck.Bk;

            /* Interpolateur temporel */
            //FIXME: antérieur à l'envoi au processeur pour éviter d'envoyer
            //tout X, trop volumineux
            Bxk = Bk(x);

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc; mp_handin;";
            
            /* Envoi des paramètres propres */
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(xos_coefs), vpack(xsd_coefs), vpack(weight), vpack(yk);} else if (mp_rank==k) {vunpack, mp_recv(0), xos_coefs; vunpack, mp_recv(0), xsd_coefs; vunpack, mp_recv(0), weight; vunpack, mp_recv(0), yk;} mp_handin;";
            
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(Bxk);} else if (mp_rank==k) {vunpack, mp_recv(0), Bxk;} mp_handin;";
        }

        /* Paramètres propres du rang 0 */
        xos_coefs = Ak_0.xos_coefs;
        xsd_coefs = Ak_0.xsd_coefs;
        weight = Ak_0.Yk.Wk.a;
        yk = Ak_0.Yk.pixels;

        /* Interpolateur temporel */
        //FIXME: antérieur à l'envoi au processeur pour éviter d'envoyer
        //tout X, trop volumineux
        Bxk = Bk_0(x);
        
        mp_exec, "if (mp_rank<nproc) {Cxk = array(double, nv); __projector, Cxk, Bxk, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, 0; Cxmyk = Cxk-yk; WCxmyk = 0.5*sum(Cxmyk*weight*Cxmyk);} else {WCxmyk=0.0;} mp_handin;";
        
        mp_exec, "if (mp_rank<nproc) {WCxk = weight*Cxk; WAxkt = array(double, nx, ny); __projector, WAxkt, WCxk, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, 1;} else {WAxkt=0.0;} mp_handin;";

        /* Calcul du gradient */
        //FIXME: le transposé de chaque projecteur a été appliqué en parallèle
        //et renvoyé au rang 0. Il reste à appliquer le transposé de
        //l'interpolateur temporel à chacune de ces rétroprojections et à
        //accumuler les résultats dans le gradient global.

        //FIXME: cette façon de faire permet à chaque job parallèle de ne
        //traiter que les parties de l'objet dont il a besoin (localité
        //temporelle des projections due à la localité de l'interpolation), ce
        //qui évite d'avoir à transférer l'intégralité de l'objet à chaque job
        //parallèle (le rapatriement serait trop coûteux en transfert mémoire)
        
        gradient = array(double, nx, ny, nt);

        /* Pour le rang 0 en premier */
        gradient(..,j_0) += w_0(-,-,)*WAxkt(..,-:1:p_0);
        
        for (k=1; k<nproc; ++k) {
            /* Interpolator */
            Bk = h_get(Cops,Ck_list(k_0+k)).Bk;
            p = Bk.wt.p;
            j = Bk.wt.j + 1;
            w = Bk.wt.w;

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc; mp_handin;";
            
            mp_exec, "if (!mp_rank) {vunpack, mp_recv(k), WAxkt;} else if (mp_rank==k) {mp_send, 0, vpack(WAxkt);} mp_handin;";

            gradient(..,j) += w(-,-,)*WAxkt(..,-:1:p);
        }
        
        mp_exec, "residuals = mp_handin(WCxmyk); mp_handin;";
        /* Residuals */
        fx += residuals;
        /* Gradient */
        gx += gradient(*);
    }
    
    return fx;
}

func trtt_2D_mpy_dyn_eval_sparse(x, &gx, Cops)
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    
    /* Extraction des paramètres communs */
    C1 = h_get(Cops,Ck_list(1));
    A1 = C1.Ak;
    B1 = C1.Bk;
    X = A1.X;
    nv = A1.Yk.nv;
    v_scl = A1.Yk.v_scl;

    /* Extraction de la fonction d'interpolation temporelle */
    funcinterp = name_of_symlink(B1.interpolator);

    /* Envoi des paramètres communs */
    mp_exec, "mp_handout, funcinterp; mp_handin;";
    /* Récupération par tous les rangs de la fonction de projection */
    mp_exec, "__interpolator = symlink_to_name(funcinterp); mp_handin;";

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
        Ak_0 = Ck_0.Ak;
        Bk_0 = Ck_0.Bk;
        p_0 = Bk_0.wt.p;
        j_0 = Bk_0.wt.j + 1;
        w_0 = Bk_0.wt.w;

        for (k=1; k<nproc; ++k) {
            Ck = h_get(Cops,Ck_list(k_0+k));
            /* Projector */
            Ak = Ck.Ak;
            xos_coefs = Ak.xos_coefs;
            xsd_coefs = Ak.xsd_coefs;
            weight = Ak.Yk.Wk.a;
            yk = Ak.Yk.pixels;
            MCk_coefs = Ak.MCk_coefs;
            MCk_row_dimlist = Ak.MCk_row_dimlist;
            MCk_row_indices = Ak.MCk_row_indices;
            MCk_col_dimlist = Ak.MCk_col_dimlist;
            MCk_col_indices = Ak.MCk_col_indices;
            /* Interpolator */
            Bk = Ck.Bk;

            /* Interpolateur temporel */
            //FIXME: antérieur à l'envoi au processeur pour éviter d'envoyer
            //tout X, trop volumineux
            Bxk = Bk(x);

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc; mp_handin;";
            
            /* Envoi des paramètres propres */
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(MCk_coefs), vpack(MCk_row_dimlist), vpack(MCk_row_indices), vpack(MCk_col_dimlist), vpack(MCk_col_indices), vpack(weight), vpack(yk);} else if (mp_rank==k) {vunpack, mp_recv(0), MCk_coefs; vunpack, mp_recv(0), MCk_row_dimlist; vunpack, mp_recv(0), MCk_row_indices; vunpack, mp_recv(0), MCk_col_dimlist; vunpack, mp_recv(0), MCk_col_indices; vunpack, mp_recv(0), weight; vunpack, mp_recv(0), yk;} mp_handin;";
            
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(Bxk);} else if (mp_rank==k) {vunpack, mp_recv(0), Bxk;} mp_handin;";
        }

        /* Paramètres propres du rang 0 */
        weight = Ak_0.Yk.Wk.a;
        yk = Ak_0.Yk.pixels;
        MCk_coefs = Ck_0.MCk_coefs;
        MCk_row_dimlist = Ck_0.MCk_row_dimlist;
        MCk_row_indices = Ck_0.MCk_row_indices;
        MCk_col_dimlist = Ck_0.MCk_col_dimlist;
        MCk_col_indices = Ck_0.MCk_col_indices;

        /* Interpolateur temporel */
        //FIXME: antérieur à l'envoi au processeur pour éviter d'envoyer
        //tout X, trop volumineux
        Bxk = Bk_0(x);
        
        mp_exec, "if (mp_rank<nproc) {Cxk = trtt_sparse_mult(Bxk, MCk_coefs, MCk_row_dimlist, MCk_row_indices, MCk_col_dimlist, MCk_col_indices, 0n); Cxmyk = Cxk-yk; WCxmyk = 0.5*sum(Cxmyk*weight*Cxmyk);} else {WCxmyk=0.0;} mp_handin;";
        
        mp_exec, "if (mp_rank<nproc) {WCxk = weight*Cxk; WAxkt = trtt_sparse_mult(WCxk, MCk_coefs, MCk_row_dimlist, MCk_row_indices, MCk_col_dimlist, MCk_col_indices, 1n);} else {WAxkt=0.0;} mp_handin;";

        /* Calcul du gradient */
        //FIXME: le transposé de chaque projecteur a été appliqué en parallèle
        //et renvoyé au rang 0. Il reste à appliquer le transposé de
        //l'interpolateur temporel à chacune de ces rétroprojections et à
        //accumuler les résultats dans le gradient global.

        //FIXME: cette façon de faire permet à chaque job parallèle de ne
        //traiter que les parties de l'objet dont il a besoin (localité
        //temporelle des projections due à la localité de l'interpolation), ce
        //qui évite d'avoir à transférer l'intégralité de l'objet à chaque job
        //parallèle (le rapatriement serait trop coûteux en transfert mémoire)
        
        gradient = array(double, nx, ny, nt);

        /* Pour le rang 0 en premier */
        gradient(..,j_0) += w_0(-,-,)*WAxkt(..,-:1:p_0);
        
        for (k=1; k<nproc; ++k) {
            /* Interpolator */
            Bk = h_get(Cops,Ck_list(k_0+k)).Bk;
            p = Bk.wt.p;
            j = Bk.wt.j + 1;
            w = Bk.wt.w;

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc; mp_handin;";
            
            mp_exec, "if (!mp_rank) {vunpack, mp_recv(k), WAxkt;} else if (mp_rank==k) {mp_send, 0, vpack(WAxkt);} mp_handin;";

            gradient(..,j) += w(-,-,)*WAxkt(..,-:1:p);
        }
        
        mp_exec, "residuals = mp_handin(WCxmyk); mp_handin;";
        /* Residuals */
        fx += residuals;
        /* Gradient */
        gx += gradient(*);
    }
    
    return fx;
}

func trtt_2D_extended_object_operator(x, job)
/* DOCUMENT trtt_2D_extended_object, x, job

   Operator for the bounds extension of the 2D + time object X. The bound conditions
   are:
   - Spatial dimensions : each temporal frame has zero values outside its
     support.
   - Temporal dimension : a periodicity condition is applied => the first frame
     follows the last.
   JOB defines if the direct or transpose of the operator is applied.
*/
{
    dims = dimsof(x);
    nx = dims(2);
    ny = dims(3);
    nt = dims(4);

    if (dims(1)==3) {
        nx = dims(2);
        ny = dims(3);
        nt = dims(4);
        if (!job) {
            out = array(double, nx+2, ny+2, nt+1);
            out(2:nx+1,2:ny+1,1:nt) = x;
            out(2:nx+1,2:ny+1,nt+1)= x(..,1);
        } else if (job==1) {
            out = x(2:nx-1,2:ny-1,1:nt-1);
            out(..,1)+=x(2:nx-1,2:ny-1,nt);
        } else {
            error, "unsupported JOB";
        }
    } else if (dims(1) == 4) {
        nx = dims(2);
        ny = dims(3);
        nz = dims(4);
        nt = dims(5);
        if (!job) {
            out = array(double, nx+2, ny+2, nz, nt+1);
            out(2:nx+1,2:ny+1,1,1:nt) = x;
            out(2:nx+1,2:ny+1,1,nt+1)= x(..,1,1);
        } else if (job==1) {
            out = x(2:nx-1,2:ny-1,1,2:nt-1);
            out(..,1,1)+=x(2:nx-1,2:ny-1,1,nt);
        } else {
            error, "unsupported JOB";
        }
    }
    return out;
}
trtt_extend = trtt_2D_extended_object_operator; // Define an easier alias

/* USER FUNCTIONS ============================================================ */

func trtt_2D_optim_cost_function_quadratic_mpy_opky(cost_mat, x, &gx)
/* DOCUMENT trtt_optim_cost_function_quadratic_mpy_opky, x, &gx, cost_mat
   
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
    nt = dims(0);
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
        fx = trtt_2D_mpy_dyn_eval(x, gx, Cops);
    } else {
        fx = trtt_2D_mpy_dyn_eval_sparse(x, gx, Cops);
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
   
    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        flag_separable = regulTV.flag_separable;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                }     
            } else if (dims(1) == 4) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wz, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                }
            }

        } else {
            /* Adding to the global cost function */
            z = trtt_extend(x,0n);
            if (!flag_separable) {
                /* Global spatio-temporal regularization */
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = trtt_extend(gx_plus,1n);
            } else {
                /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = trtt_extend(gx_plus,1n);
            }
        }
        
        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    }
    
    else if (!is_void(regul)) {

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
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

func trtt_2D_optim_cost_function_quadratic_opky(cost_mat, x, &gx)
/* DOCUMENT trtt_optim_cost_function_quadratic_opky, cost_mat, x, &gx
   
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
    nt = dims(0);
    /* Projector */
    Cops = cost_mat.Cops;
   
    /* True data transposition */
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

    Ck_list = h_get(Cops,"Ck_list");
    ndata = numberof(Ck_list);
    fx = 0.0;
    gx = 0.0;
    for (k=1; k<=ndata; ++k) {
        /* Get the projector */
        Ck = h_get(Cops, Ck_list(k));
        Ak = Ck.Ak;
        Bk = Ck.Bk;
        Yk = Ak.Yk;
        /* Get corresponding data */
        yk = Yk.pixels;
        /* Get noise covariance matrix */
        Wk = Yk.Wk;
        
        if (!exist_sparse_coefs) {
            Cxk = Ck(x);
        } else {
            MCk = Ak.MCk;
            Cxk = MCk((Bk(x))(*));
        }
  
        Cxmyk = Cxk-yk;
        WCxmyk = Wk(Cxmyk);
        fx += 0.5*sum(Cxmyk*WCxmyk);

        if (!exist_sparse_coefs) {
            gx += Ck(Wk(Cxk),1)(*);
        } else {
            MCk = Ak.MCk;
            gx += Bk(MCk(Wk(Cxk),1),1)(*);
        }
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
    
    /* Adding a constraint function for regularization */ 
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        flag_separable = regulTV.flag_separable;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                }     
            } else if (dims(1) == 4) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wz, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                }
            }

        } else {
            /* Adding to the global cost function */
            z = trtt_extend(x,0n);
            if (!flag_separable) {
                /* Global spatio-temporal regularization */
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = trtt_extend(gx_plus,1n);
            } else {
                /* Separable spatio-temporal regularization */
                if (numberof(threshold)==1) {
                    thr_t = threshold;
                    thr_s = threshold;
                } else {
                    thr_t = threshold(2);
                    thr_s = threshold(1);
                }     
                fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                gx_plus = trtt_extend(gx_plus,1n);
            }
        }
        
        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    }
    
    else if (!is_void(regul)) {

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
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
trtt_cost_quadratic_opky = trtt_2D_optim_cost_function_quadratic_opky; // Define an easier alias

func trtt_2D_optim_cost_function_quadratic_mpy(x, &gx, cost_mat)
/* DOCUMENT trtt_optim_cost_function_quadratic_mpy, x, &gx, cost_mat
   
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
    nt = dims(0);
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
        fx = trtt_2D_mpy_dyn_eval(x, gx, Cops);
    } else {
        fx = trtt_2D_mpy_dyn_eval_sparse(x, gx, Cops);
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
   
    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        flag_separable = regulTV.flag_separable;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                }     
            } else if (dims(1) == 4) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wz, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                }
            }

        } else {
            /* Adding to the global cost function */
            z = trtt_extend(x,0n);
            if (!flag_separable) {
                /* Global spatio-temporal regularization */
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = trtt_extend(gx_plus,1n);
            } else {
                /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = trtt_extend(gx_plus,1n);
            }
        }
        
        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    }
    
    else if (!is_void(regul)) {

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
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

func trtt_2D_optim_cost_function_quadratic(x, &gx, cost_mat)
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
    nt = dims(0);
    /* Projector */
    Cops = cost_mat.Cops;
   
    /* True data transposition */
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

    Ck_list = h_get(Cops,"Ck_list");
    ndata = numberof(Ck_list);
    fx = 0.0;
    gx = 0.0;
    for (k=1; k<=ndata; ++k) {
        /* Get the projector */
        Ck = h_get(Cops, Ck_list(k));
        Ak = Ck.Ak;
        Bk = Ck.Bk;
        Yk = Ak.Yk;
        /* Get corresponding data */
        yk = Yk.pixels;
        /* Get noise covariance matrix */
        Wk = Yk.Wk;
        
        if (!exist_sparse_coefs) {
            Cxk = Ck(x);
        } else {
            MCk = Ak.MCk;
            Cxk = MCk((Bk(x))(*));
        }
  
        Cxmyk = Cxk-yk;
        WCxmyk = Wk(Cxmyk);
        fx += 0.5*sum(Cxmyk*WCxmyk);

        if (!exist_sparse_coefs) {
            gx += Ck(Wk(Cxk),1)(*);
        } else {
            MCk = Ak.MCk;
            gx += Bk(MCk(Wk(Cxk),1),1)(*);
        }
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
    
    /* Adding a constraint function for regularization */ 
    if (!is_void(regulTV)) {

        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        flag_separable = regulTV.flag_separable;
        local gx_plus;

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wt, 1);
                }     
            } else if (dims(1) == 4) {
                /* Adding to the global cost function */
                z = trtt_extend(spl_interp_apply(x, wx, wy, wz, wt),0n);
                if (!flag_separable) {
                    /* Global spatio-temporal regularization */
                    fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                } else {
                    /* Separable spatio-temporal regularization */
                    if (numberof(threshold)==1) {
                        thr_t = threshold;
                        thr_s = threshold;
                    } else {
                        thr_t = threshold(2);
                        thr_s = threshold(1);
                    }     
                    fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                    gx_plus = spl_interp_apply(trtt_extend(gx_plus,1n), wx, wy, wz, wt, 1);
                }
            }

        } else {
            /* Adding to the global cost function */
            z = trtt_extend(x,0n);
            if (!flag_separable) {
                /* Global spatio-temporal regularization */
                fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
                gx_plus = trtt_extend(gx_plus,1n);
            } else {
                /* Separable spatio-temporal regularization */
                if (numberof(threshold)==1) {
                    thr_t = threshold;
                    thr_s = threshold;
                } else {
                    thr_t = threshold(2);
                    thr_s = threshold(1);
                }     
                fx_plus = rgl_mixed_ndpt(weight(1), thr_s, weight(2), thr_t, z, gx_plus, 0n);
                gx_plus = trtt_extend(gx_plus,1n);
            }
        }
        
        /* Adding to the global cost function */
        gx += gx_plus;
        fx += fx_plus;
    }

    // else if (!is_void(regulTVsep)) { //FIXME: TEST TOTAL VARIATION SEPARABLE

    //     weight = regulTVsep.weight;
    //     threshold = regulTVsep.threshold;
    //     options = regulTVsep.options;
    //     local gx_plus;
    //     local gx_plus2;

    //     if (type == TRTT_SPLINE_DRIVEN) {        
    //         /* For regularization we have to work in samples space, and not in
    //          * spline coefficients space */
    //         // if (dims(1) == 3) {              
    //             /* Adding to the global cost function */
    //             z = spl_interp_apply(x, wx, wy, wt);
    //             /* spatial part */
    //             fx_plus = rgl_totvar(z, gx_plus, weight=[weight(1),weight(2),0.0], threshold=threshold, options=options);                               
    //             gx_plus = spl_interp_apply(gx_plus, wx, wy, wt, 1);
    //             /* temporal part */
    //             fx_plus2 = rgl_totvar(z, gx_plus2, weight=[0.0,0.0,weight(3)], threshold=threshold, options=options);                               
    //             gx_plus2 = spl_interp_apply(gx_plus2, wx, wy, wt, 1);
    //             fxplus += fxplus2;
    //             gxplus += gxplus2;          
    //         // } else if (dims(1) == 4) {
    //         //     /* Adding to the global cost function */
    //         //     z = spl_interp_apply(x, wx, wy, wz, wt);
    //         //     fx_plus = rgl_totvar(z, gx_plus, weight=weight, threshold=threshold, options=options);
    //         //     gx_plus = spl_interp_apply(gx_plus, wx, wy, wz, wt, 1);
    //         // }

    //     // } else {
    //     //     /* Adding to the global cost function */
    //     //     fx_plus = rgl_totvar(x, gx_plus, weight=weight, threshold=threshold, options=options);
    //     // }

    //     /* Adding to the global cost function */
    //     gx += gx_plus;
    //     fx += fx_plus;
    // }
  
    else if (!is_void(regul)) {

        if (type == TRTT_SPLINE_DRIVEN) {        
            /* For regularization we have to work in samples space, and not in
             * spline coefficients space */
            if (dims(1) == 3) {
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
trtt_cost_quadratic = trtt_2D_optim_cost_function_quadratic; // Define an easier alias

func trtt_2D_optim_define_cost_mat(Cops, dweights, use_sparse_coefs=, matrix=, regul=, regulTV=, viewer=, win_viewer=, win_viewer2=, cmin=, cmax=, roi=, xref=)
/* DOCUMENT CT_system_OptimPack_define_cost_mat, Cops, dweights, use_sparse_coefs=, regul=, regulTV=, viewer=, win_viewer=, win_viewer2=, cmin=, cmax=, roi=, xref=
   
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
    Ak1 = C1.Ak;
    X = Ak1.X;
    nx = X.nx;
    ny = X.ny;
    nt = X.nt;
    nv = Ak1.Yk.nv;

    cost_mat = h_new(
        Cops=Cops,
        regul=regul,
        regulTV=regulTV,
        viewer=viewer,
        win_viewer=win_viewer,
        win_viewer2=win_viewer2,
        cmin=cmin,
        cmax=cmax,
        nt=nt,
        dyn=1n);

    /* Calculate true data transposition */
    WCyT = 0.0;
    if (use_sparse_coefs & !Cops.exist_sparse_coefs) {
        _trtt_2D_sparse_coefs, Cops, dyn=1n;
    }
    for (k=1; k<=ndata; ++k) {
        Ck = h_get(Cops, Ck_list(k));
        Ak = Ck.Ak;
        if (Cops.exist_sparse_coefs) {
            MCk = sparse_matrix(Ak.MCk_coefs, Ak.MCk_row_dimlist, Ak.MCk_row_indices, Ak.MCk_col_dimlist, Ak.MCk_col_indices);
            h_set, Ak, MCk=MCk;
        }
        Yk = Ak.Yk;
        yk = Yk.pixels;
        weight = Yk.weight;
        if (!is_void(weight)) {
            Wk = linop_new(LINOP_DIAGONAL, dweights(,k)*weight);
        } else {
            Wk = linop_new(LINOP_DIAGONAL, dweights(,k)*array(1.0,numberof(yk)));
        }    
        WCyT += Ck(Wk(yk),1);
        h_set, Yk, Wk = Wk;
    }
    /* Set data transposition in COST_MAT */
    h_set, cost_mat, WCyT = WCyT;
    
    /* Set interpolation coefficients for SPLINE DRIVEN only */
    type = Cops.type;
    
    if (type == TRTT_SPLINE_DRIVEN) {
        h_set, cost_mat, wx = X.wx;
        h_set, cost_mat, wy = X.wy;
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

func trtt_2D_optim_simu_launcher(Cops, cost, x=, method=, mem=, frtol=, fatol=, maxeval=, dweights=, get_costs=, metrics=, roi=, xref=, use_sparse_coefs=, matrix=, regul=, regulTV=, xmin=, xmax=, maxiter=, viewer=, win_viewer=, win_viewer2=, cmin=, cmax=, verbose=, store=)        
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
    X = h_get(Cops,Ck_list(1)).Ak.X;
    
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
        nv = h_get(Cops,Ck_list(1)).Ak.Yk.nv;
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
        if (cost==trtt_cost_quadratic || trtt_cost_quadratic_mpy) {
            printer = _trtt_2D_metrics;
        } else if (cost==trtt_cost_quadratic_opky || cost==trtt_cost_quadratic_mpy_opky) {
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
        if (cost==trtt_cost_quadratic || cost==trtt_cost_quadratic_mpy) {
            printer=std_viewer;
        } else if (cost==trtt_cost_quadratic_opky || cost==trtt_cost_quadratic_mpy_opky) {
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
    if (cost==trtt_cost_quadratic || cost==trtt_cost_quadratic_mpy) {
        mem_copy, &x, op_mnb(cost, x, extra=cost_mat, xmin=xmin, xmax=xmax, method=method, mem=mem, maxiter=maxiter, maxeval=maxeval, frtol=frtol, fatol=fatol, viewer=printer, verb=verbose);
    } else if (cost==trtt_cost_quadratic_opky) { 
        /* FIXME: ADAPTATION POUR NOUVEL OPTIMPACK OPKY */
        minimizer=closure("trtt_cost_quadratic_opky",cost_mat);
        mem_copy, &x, opk_minimize(minimizer, x, mem=mem, lower=xmin, upper=xmax, maxiter=maxiter, maxeval=maxeval, printer=printer, verb=verbose);
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
        mean_noise_std += Ck.Ak.Yk.mean_var_noise;
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

// func trtt_2D_optim_standard_viewer(x, cost_mat)
// /* DOCUMENT trtt_optim_standard_viewer, x, cost_mat
   
//    Viewer for the call of the optimization algorithm for the reconstruction of a
//    phantom. Plot the current estimate X.

// */
// {
//     nwin = cost_mat.win_viewer;
//     cmin = cost_mat.cmin;
//     cmax = cost_mat.cmax;

//     if (is_void(cmin)) cmin = min(x(*));
//     if (is_void(cmax)) cmax = max(x(*));

//     if (is_void(nwin)) nwin=0;
        
//     mywindow, nwin; fma; img_plot, x(..,1), cbar=1, cmin=cmin, cmax=cmax, vert=1, format="%.3f";

//     pause, 1;
// }
// std_viewer = trtt_2D_optim_standard_viewer; // Define an easier alias

// func trtt_2D_optim_low_contrast_viewer(x, cost_mat)
// /* DOCUMENT trtt_optim_low_constrast_viewer, x, cost_mat
   
//    Viewer for the call of the optimization algorithm for the reconstruction of a
//    low contrast phantom (the image X will be histogram equalized), for example
//    the Shepp Logan phantom. Plot the current estimate X.
// */
// {
//     nwin = cost_mat.win_viewer;
//     cmin = cost_mat.cmin;
//     cmax = cost_mat.cmax;

//     if (is_void(cmin)) cmin = min(x);
//     if (is_void(cmax)) cmax = max(x);

//     if (is_void(nwin)) nwin=0;
        
//     trtt_2D_egalize_and_display_low_contrast, x(..,1), cmin=cmin, cmax=cmax, nwin=nwin;

//     pause, 1;
// }
// low_ctr_viewer = trtt_2D_optim_low_contrast_viewer; // Define an easier alias

// func trtt_2D_egalize_and_display_low_contrast(z, first=, scale=, top=, cmin=, cmax=, nticks=, nwin=)
// /* DOCUMENT trtt_egalize_and_display_low_contrast, z, top=, cmin=, cmax=, nticks=, nwin=
   
//    This function is based on the HISTEQ_SCALE function which applies an
//    egalization of the histogram of Z. Then the resulting image is plotted with
//    the appropriate colorbar involving the "egalized" levels.  NTICKS= : Number
//    of plotted labels on the colorbar.  NWIN= : Number of the display
//    window. Default is 0.
   
//    SEE ALSO: histeq_scale, img_plot.
// */
// {
//     if (is_void(nticks)) nticks = 11;
//     if (is_void(nwin)) nwin=0;
//     if (is_void(cmin)) cmin = min(z(*));
//     if (is_void(cmax)) cmax= max(z(*));

        
//     if (is_void(top)) top= bytscl([0.,1.])(2);  /* palette size - 1 */
//     top= long(top);
//     if (top<0 | top>255) error, "top value out of range 0-255";
//     y= z(*);
//     if (!is_void(cmin)) y= y(where(y>=cmin));
//     if (!is_void(cmax)) y= y(where(y<=cmax));
//     y= y(sort(y));
//     x= span(0.,1., numberof(y));
//     xp= span(0.,1., top+2);
//     bins= interp(y, x, xp);
//     list= where(bins(dif)<=0.0);
//     if (numberof(list)) {
//         /* some value (or values) of z are repeated many times --
//            try to handle this by adding a small slope to the sorted y */
//         dy= y(0)-y(1);
//         if (!dy) dy= 1.0;
//         for (eps=1.e-10 ; eps<1000.1 ; eps*=10.) {
//             bins= interp(y+eps*dy*x, x, xp);
//             list= where(bins(dif)<=0.0);
//             if (!numberof(list)) break;
//         }
//         if (eps>1000.) error, "impossible error??";
//     }

//     // Apply the contrast egalization
//     zs = char(max(min(digitize(z,bins)-2,top),0));

//     // Define the labels
//     levstep = int((top+2)/nticks)+1;
//     labels = bins(1:0:levstep);
       
//     // Plot the egalized image with the appropriate colorbar
//     mywindow, nwin; fma; img_plot, zs, first=first, scale=scale, cbar=1, cmin=bins(1), cmax=bins(0), vert=1, nticks=nticks, format="%.3f", labels=labels;
// }
// trtt_pli_low_ctr = trtt_2D_egalize_and_display_low_contrast; // define an easier alias

// func trtt_2D_mpy_dyn_eval_old(x, &gx, Cops)
// {
//     Ck_list = Cops.Ck_list;
//     ndata = numberof(Ck_list);
    
//     /* Extraction des paramètres communs */
//     C1 = h_get(Cops,Ck_list(1));
//     A1 = C1.Ak;
//     B1 = C1.Bk;
//     X = A1.X;
//     nx = X.nx; ny = X.ny; nt = X.nt;
//     s_scl = X.s_scl;
//     s_deg = X.s_deg;
//     nv = A1.Yk.nv;
//     v_scl = A1.Yk.v_scl;
//     local footprint; eq_nocopy, footprint, X.footprint;
//     local cum_footprint; eq_nocopy, cum_footprint, X.cum_footprint;
//     size_footprint = X.size_footprint;
//     step_footprint = X.step_footprint;
//     local coord_footprint; eq_nocopy, coord_footprint, X.coord_footprint;

//     /* Extraction de la fonction de projection */
//     funcproj = name_of_symlink(A1.projector);
//     /* Extraction de la fonction d'interpolation temporelle */
//     funcinterp = name_of_symlink(B1.interpolator);

//     /* Envoi des paramètres communs */
//     mp_exec, "mp_handout, x, nx, ny, nt, s_scl, s_deg, nv, v_scl, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, funcproj, funcinterp;";
//     /* Récupération par tous les rangs de la fonction de projection */
//     mp_exec, "__projector = symlink_to_name(funcproj);";
//     mp_exec, "__interpolator = symlink_to_name(funcinterp);";

//     /* Calcul du nombre d'occurrences */
//     Nmult = ndata/mp_size;
//     Nrest = ndata%mp_size;

//     /* Initialisation des résidus */
//     fx = 0.0;

//     /* Calcul des occurrences */
//     for (i=0; i<=Nmult; ++i) {

//         if (i<Nmult) nproc = mp_size;
//         else if (Nrest > 0) nproc = Nrest;
//         else break;
        
//         k_0 = i*mp_size + 1;
//         Ck_0 = h_get(Cops,Ck_list(k_0));
//         Ak_0 = Ck_0.Ak;
//         Bk_0 = Ck_0.Bk;

//         for (k=1; k<nproc; ++k) {
//             Ck = h_get(Cops,Ck_list(k_0+k));
//             /* Projector */
//             Ak = Ck.Ak;
//             xos_coefs = Ak.xos_coefs;
//             xsd_coefs = Ak.xsd_coefs;
//             weight = Ak.Yk.Wk.a;
//             yk = Ak.Yk.pixels;
//             /* Interpolator */
//             Bk = Ck.Bk;
//             p = Bk.wt.p;
//             j = Bk.wt.j + 1;
//             w = Bk.wt.w;

//             /* Envoi du numéro de rang concerné + */
//             mp_exec, "mp_handout, k, nproc";
            
//             /* Envoi des paramètres propres */
//             mp_exec, "if (!mp_rank) {mp_send, k, vpack(xos_coefs), vpack(xsd_coefs), vpack(weight), vpack(yk);} else if (mp_rank==k) {vunpack, mp_recv(0), xos_coefs; vunpack, mp_recv(0), xsd_coefs; vunpack, mp_recv(0), weight; vunpack, mp_recv(0), yk;} ";
            
//             mp_exec, "if (!mp_rank) {mp_send, k, p, vpack(j), vpack(w);} else if (mp_rank==k) {p = mp_recv(0); vunpack, mp_recv(0), j; vunpack, mp_recv(0), w;} ";
//         }

//         /* Paramètres propres du rang 0 */
//         xos_coefs = Ak_0.xos_coefs;
//         xsd_coefs = Ak_0.xsd_coefs;
//         weight = Ak_0.Yk.Wk.a;
//         yk = Ak_0.Yk.pixels;
//         p = Bk_0.wt.p;
//         j = Bk_0.wt.j + 1;
//         w = Bk_0.wt.w;

//         mp_exec, "if (mp_rank<nproc) {Bxk = array(double, nx, ny); __interpolator, Bxk, x, p, j, w, 0; Cxk = array(double, nv); __projector, Cxk, Bxk, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, 0; Cxmyk = Cxk-yk; WCxmyk = 0.5*sum(Cxmyk*weight*Cxmyk);} else {WCxmyk=0.0;}";
        
//         mp_exec, "if (mp_rank<nproc) {WCxk = weight*Cxk; WAxkt = array(double, nx, ny); __projector, WAxkt, WCxk, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, 1; WCxkt = array(double, nx, ny, nt); __interpolator, WCxkt, WAxkt, p, j, w, 1;} else {WCxkt = 0.0;}";
        
//         mp_exec, "residuals = mp_handin(WCxmyk); gradient = mp_handin(WCxkt(*));";
//         /* Residuals */
//         fx += residuals;
//         /* Gradient */
//         gx += gradient;
//     }
    
//     return fx;
// }


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
