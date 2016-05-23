/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_metrics.i --
 *
 * TRTT 2D evaluation metrics package.
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

extern trtt_3D_metrics_info
/* DOCUMENT trtt_3D_metrics_info

   When a reconstruction is performed with the function
   trtt_3D_optim_simu_launcher(), we must specify, with the parameter METRICS, a
   series of evaluation metrics to calculate for the quantification of
   results. METRICS is an table of booleans each element of which answers
   "yes/no" to the question : "Do we want to calculate this metric ?". Here is
   the table showing which element correponds to what metric in a predefined
   order. We must give the whole array with each element in this order.
 
    _______
   |       |
      0/1     =>    1 : residuals at the end of reconstruction.

      0/1     =>    2 : residuals at each iteration of reconstruction.

      0/1     =>    3 : regularization term at the end of reconstruction.

      0/1     =>    4 : regularization term at each iteration of reconstruction.

      0/1     =>    5 : RMS error at the end of reconstruction (need XREF).

      0/1     =>    6 : RMS error at each iteration of reconstruction (need XREF).

   |_______|

   Example : In a given reconstruction we want to calculate the residuals at
   each iteration, the regularization term at the end and the RMS error at the
   end. The parameter METRICS in trtt_3D_optim_simu_launcher() will be as
   follows :

                                [0, 1, 1, 0, 1, 0]

   SEE ALSO: trtt_3D_optim_simu_launcher.
   
 */

/* GLOBALS =================================================================== */

TRTT_3D_METRICS_NUMBER = 6;
/* DOCUMENT number of available metrics.

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

/* FUNCTIONS ================================================================= */

func _trtt_3D_metrics(x, cost_mat)
/* DOCUMENT _trtt_3D_metrics, x, cost_mat
   
   Calculate several metrics when called at each iteration of a reconstruction
   process. It is used in OP_MNB as a "markup" function.

   SEE ALSO: op_mnb, trtt_3D_metrics_info.
 */
{
    extern TRTT_3D_RESIDUALS_ITER;
    extern TRTT_3D_REGULTERM_ITER;
    extern TRTT_3D_RMSE_ITER;

    extern rec_results;

    viewer = cost_mat.viewer;
    
    if (TRTT_3D_RESIDUALS_ITER) {
        local residuals_iter; eq_nocopy, residuals_iter, rec_results.residuals_iter;
        grow, residuals_iter, _trtt_3D_residuals(x, cost_mat.Cops, dyn=cost_mat.dyn);
        h_set, rec_results, residuals_iter=residuals_iter;
    }
    if (TRTT_3D_REGULTERM_ITER) {
        local regulterm_iter; eq_nocopy, regulterm_iter, rec_results.regulterm_iter;
        grow, regulterm_iter, _trtt_3D_regulterm(x, regul=cost_mat.regul, regulTV=cost_mat.regulTV);
        h_set, rec_results, regulterm_iter=regulterm_iter;
    }
    if (TRTT_3D_RMSE_ITER) {
        Cops = cost_mat.Cops;
        local rmse_iter; eq_nocopy, rmse_iter, rec_results.rmse_iter;
        if (Cops.type == TRTT_SPLINE_DRIVEN) {
            if (dims(1) == 3)
                x = spl_interp_apply(x, cost_mat.wx, cost_mat.wy, cost_mat.wz);
            else
                x = spl_interp_apply(x, cost_mat.wx, cost_mat.wy, cost_mat.wz, cost_mat.wt);
        }
        grow, rmse_iter, _trtt_3D_rmse(x, cost_mat.xref, cost_mat.s_deg, cost_mat.t_deg, roi=cost_mat.roi);
        h_set, rec_results, rmse_iter=rmse_iter;
    }

    /* Apply viewer function */
    if (viewer) {
        std_viewer, x, cost_mat;
    }
}

func trtt_3D_optim_standard_viewer(x, cost_mat)
/* DOCUMENT trtt_3D_optim_standard_viewer, x, cost_mat
   
   Viewer for the call of the optimization algorithm for the reconstruction of a
   phantom. Plot the current estimate X.

*/
{
    nz = cost_mat.nz;
    kz = long(floor(nz/2.0));
    nwin = cost_mat.win_viewer;
    cmin = cost_mat.cmin;
    cmax = cost_mat.cmax;

    if (is_void(cmin)) cmin = min(x(*));
    if (is_void(cmax)) cmax = max(x(*));

    if (is_void(nwin)) nwin=0;

    if (!(cost_mat.dyn)) {
        trtt3D_plot_vox, x, 3, kz, nwin, cmin=cmin, cmax=cmax;
    } else {
        idt=int(floor(0.5*cost_mat.nt));
        trtt3D_plot_slice, x, 3, kz, idt, nwin, cmin=cmin, cmax=cmax;
    }
    cmap, "gray";

    pause, 1;
}
std_viewer = trtt_3D_optim_standard_viewer; // Define an easier alias

func _trtt_3D_metrics_for_opky(cost_mat, x, fx, gx, iter, eval, t)
/* DOCUMENT _trtt_3D_metrics_for_opky, cost_mat, x, fx, gx, iter, eval, t
   
   Calculate several metrics when called at each iteration of a reconstruction
   process. It is used in OP_MNB as a "markup" function.

   SEE ALSO: opk_minimize, trtt_3D_metrics_info.
 */
{
    extern TRTT_3D_RESIDUALS_ITER;
    extern TRTT_3D_REGULTERM_ITER;
    extern TRTT_3D_RMSE_ITER;

    extern rec_results;

    viewer = cost_mat.data.viewer;
    
    if (TRTT_3D_RESIDUALS_ITER) {
        local residuals_iter; eq_nocopy, residuals_iter, rec_results.residuals_iter;
        grow, residuals_iter, _trtt_3D_residuals(x, cost_mat.data.Cops, dyn=cost_mat.data.dyn);
        h_set, rec_results, residuals_iter=residuals_iter;
    }
    if (TRTT_3D_REGULTERM_ITER) {
        local regulterm_iter; eq_nocopy, regulterm_iter, rec_results.regulterm_iter;
        grow, regulterm_iter, _trtt_3D_regulterm(x, regul=cost_mat.data.regul, regulTV=cost_mat.data.regulTV);
        h_set, rec_results, regulterm_iter=regulterm_iter;
    }
    if (TRTT_3D_RMSE_ITER) {
        Cops = cost_mat.data.Cops;
        local rmse_iter; eq_nocopy, rmse_iter, rec_results.rmse_iter;
        if (Cops.type == TRTT_SPLINE_DRIVEN) {
            if (dims(1) == 3)
                x = spl_interp_apply(x, cost_mat.data.wx, cost_mat.data.wy, cost_mat.data.wz);
            else
                x = spl_interp_apply(x, cost_mat.data.wx, cost_mat.data.wy, cost_mat.data.wz, cost_mat.data.wt);
        }
        grow, rmse_iter, _trtt_3D_rmse(x, cost_mat.data.xref, cost_mat.data.s_deg, cost_mat.data.t_deg, roi=cost_mat.data.roi);
        h_set, rec_results, rmse_iter=rmse_iter;
    }

    /* Apply viewer function */
    if (viewer) {
        std_viewer_opky, cost_mat, x, fx, gx, iter, eval, t;
    }
}

func trtt_3D_optim_standard_viewer_for_opky(cost_mat, x, fx, gx, iter, eval, t)
/* DOCUMENT trtt_3D_optim_standard_viewer_for_opky, cost_mat, x, fx, gx, iter, eval, t
   
   Viewer for the call of the optimization algorithm for the reconstruction of a
   phantom. Plot the current estimate X.

*/
{
    nz = cost_mat.data.nz;
    kz = long(floor(nz/2.0));
    nwin = cost_mat.data.win_viewer;
    nwin2 = cost_mat.data.win_viewer2;
    cmin = cost_mat.data.cmin;
    cmax = cost_mat.data.cmax;

    if (is_void(cmin)) cmin = min(x(*));
    if (is_void(cmax)) cmax = max(x(*));

    if (is_void(nwin)) nwin=0;

    if (!(cost_mat.data.dyn)) {
        trtt3D_plot_vox, x, 3, kz, nwin, cmin=cmin, cmax=cmax;
    } else {
        idt=int(floor(0.5*cost_mat.data.nt));
        trtt3D_plot_slice, x, 3, kz, idt, nwin, cmin=cmin, cmax=cmax;
    }
    cmap, "gray";

    if (!is_void(nwin2)) {
        window, nwin2; plg, fx, iter, marks=1, marker='\2', color="black", width=4.0, type="none";
        logxy, 0, 1;
        limits, "e", "e", fx/10.0, "e";
    }

    pause, 1;
}
std_viewer_opky = trtt_3D_optim_standard_viewer_for_opky; // Define an easier alias

func _trtt_3D_rmse(x, xref, s_deg, t_deg, roi=)
/* DOCUMENT _trtt_3D_rmse, x, xref, s_deg, t_deg, roi=
   
   Calculate the root mean square error between X and XREF, which are
   spatio-temporal 3D+t images. The error can be calculated on a spatial region
   of interest (ROI). Thus ROI is an 4 element array corresponding to the bounds
   (indexes) of the ROI : [id_min,idx_max,idy_min,idy_max].
 */
{
    dims = dimsof(x);
    dimsref = dimsof(xref);
    nx = dims(2); ny = dims(3); nz = dims(4);
    nxref = dimsref(2); nyref = dimsref(3); nzref = dimsref(4);

    if (dims(1) == 3) {
        w = spl_spline_interpolator(indgen(0:nxref-1)*(nx-1.)/(nxref-1.), nx, s_deg, indgen(0:nyref-1)*(ny-1.)/(nyref-1.), ny, s_deg, indgen(0:nzref-1)*(nz-1.)/(nzref-1.), nz, s_deg);
    } else if (dims(1) == 4) {
        nt = dims(5); ntref = dimsref(5);
        w = spl_spline_interpolator(indgen(0:nxref-1)*(nx-1.)/(nxref-1.), nx, s_deg, indgen(0:nyref-1)*(ny-1.)/(nyref-1.), ny, s_deg, indgen(0:nzref-1)*(nz-1.)/(nzref-1.), nz, s_deg, indgen(0:ntref-1)*(nt-1.)/(ntref-1.), nt, t_deg);
    }

    if (!is_void(roi)) {
        if (dims(1) == 3) {
            return (w(spl_spline_samples_to_coefficients(x,s_deg,s_deg,s_deg))-xref)(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6))(*)(rms);
        } else if (dims(1) == 4) {
            return (w(spl_spline_samples_to_coefficients(x,s_deg,s_deg,s_deg,t_deg))-xref)(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6),..)(*)(rms);
        }
    } else {
        if (dims(1) == 3) {
            return (w(spl_spline_samples_to_coefficients(x,s_deg,s_deg,s_deg))-xref)(*)(rms);
        } else if (dims(1) == 4) {
            return (w(spl_spline_samples_to_coefficients(x,s_deg,s_deg,s_deg,t_deg))-xref)(*)(rms);
        }
    }
}  

func _trtt_3D_residuals(x, Cops, dyn=)
/* DOCUMENT _trtt_3D_residuals, x, Cops, dyn=

   Calculate the residuals of a reconstruction, given the reconstructed image X,
   the projector COPS which contains the data Y.
 */
{
    Ck_list = h_get(Cops,"Ck_list");
    ndata = numberof(Ck_list);

    residuals = 0.0;
    for (k=1; k<=ndata; ++k) {
        /* Get the projector */
        Ck = h_get(Cops, Ck_list(k));
        /* Get corresponding data */
        if (!dyn) {
            yk = Ck.Yk.pixels;
            nv = Ck.Yk.nv;
            nu = Ck.Yk.nu;
        } else {
            yk = Ck.Ak.Yk.pixels;
            nv = Ck.Ak.Yk.nv;
            nu = Ck.Ak.Yk.nu;
        }
        
        Cxk = Ck(x);
        Cxmyk = Cxk-yk;
        /* Calculate residuals */
        residuals += (1./(nu*nv))*sum(Cxmyk*Cxmyk);
    }
    residuals /= double(ndata);
    /* Return residuals */
    return residuals;
}

func _trtt_3D_regulterm(x, regul=, regulTV=)
/* DOCUMENT _trtt_3D_regulterm, x, regul, regulTV=

   Calculate the regularization term of a reconstruction, given the
   reconstructed image X, and the regularization REGUL. If not given the
   regularization term will be 0.
 */
{
    dims = dimsof(x);

    /* Adding a constraint function for regularization */
    if (!is_void(regulTV)) {
        weight = regulTV.weight;
        threshold = regulTV.threshold;
        options = regulTV.options;
        regulterm = rgl_totvar(x, weight=weight, threshold=threshold, options=options);
    } else if (!is_void(regul)) {
        // Updating the regularization term
        rgl_update, regul, x;
        // Adding to the global cost function
        regulterm = rgl_get_penalty(regul, x);
    } else {
        regulterm = 0.0;
    }
        
    /* Return regularization term */
    return regulterm;
}
