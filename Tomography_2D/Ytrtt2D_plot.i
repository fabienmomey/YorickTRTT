/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt2D_plot.i --
 *
 * TRTT 2D plot function package.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011, the MiTiV Team.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt2D_plot.i --
 *
 * TRTT 2D plot function package.
 *
 *--------------------------------------------------------------
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

/* USER FUNCTIONS ============================================================ */

func trtt_2D_get_sinogram(Cops, dyn=)
/* DOCUMENT trtt_2D_get_sinogram, Cops, dyn=

   get sinogram data stored in the projectors set COPS.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    if (!dyn) {
        nv = h_get(Cops,Ck_list(1)).Yk.nv;
    } else {
        nv = h_get(Cops,Ck_list(1)).Ak.Yk.nv;
    }
    sinogram=array(double,nv,ndata);

    for (k=1; k<=ndata; ++k) {
        if (!dyn) {
            Yk = h_get(Cops, Ck_list(k)).Yk;
        } else {
            Yk = h_get(Cops, Ck_list(k)).Ak.Yk;
        }
        sinogram(,k) = Yk.pixels;
    }
    
    return sinogram;
}
trtt_get_sino = trtt_2D_get_sinogram;

func trtt_2D_get_ref_sinogram(Htomo, sparse=)
/* DOCUMENT trtt_2D_get_ref_sinogram, Htomo

   Get projections of the reference voxel image XREF "passed through" the
   projectors Ck contained in COPS.
 */
{
    if (is_void(sparse)) sparse=0n;
    X = Htomo.X;
    Y = Htomo.Y;
    Xref = Htomo.Xref;
    Cops = Htomo.Cops;
    exist_sparse_coefs = Cops.exist_sparse_coefs;
    s_deg = X.s_deg;
    t_deg = X.t_deg;
    x = Xref.voxels;
    dims=dimsof(x);
    if (type == TRTT_SPLINE_DRIVEN) {
        if (dims(1)==2)
            c = spl_spline_samples_to_coefficients(x, s_deg, s_deg);
        else if (dims(1)==3)
            c = spl_spline_samples_to_coefficients(x, s_deg, s_deg, t_deg);
        else
            error, "Bad number of dimensions.";
    }
    Ck_list = Cops.Ck_list;
    Yk_list = Y.Yk_list;
    ndata = numberof(Ck_list);
    nv = h_get(Y,Yk_list(1)).nv;
    ya = array(double,nv,ndata);
    for (k=1; k<=ndata; ++k) {
        write, format="Projection n°%d\n", k;
        Ck = h_get(Cops,Ck_list(k));
        Yk = h_get(Y,Yk_list(k));
        if (type == TRTT_SPLINE_DRIVEN)
            if (!sparse || !exist_sparse_coefs) {
            ya(,k) = Ck(c);
            } else {
                if (is_void(Ck.MCk))
                    MCk = sparse_matrix(Ck.MCk_coefs, Ck.MCk_row_dimlist, Ck.MCk_row_indices, Ck.MCk_col_dimlist, Ck.MCk_col_indices);
                else
                    MCk = Ck.MCk;
                ya(,k) = MCk(c(*));
            }
        else
            if (!sparse || !exist_sparse_coefs) {
            ya(,k) = Ck(x);
            } else {
                if (is_void(Ck.MCk))
                    MCk = sparse_matrix(Ck.MCk_coefs, Ck.MCk_row_dimlist, Ck.MCk_row_indices, Ck.MCk_col_dimlist, Ck.MCk_col_indices);
                else
                    MCk = Ck.MCk;
                ya(,k) = MCk(x(*));
            }
    }

    return ya;
}
trtt_get_ref_sino = trtt_2D_get_ref_sinogram;

func trtt_2D_plot_tomobj_slice(X, idt, nwin, cmin=, cmax=)
/* DOCUMENT trtt_2D_plot_tomobj_slice, X, idt, nwin, cmin=, cmax=

   Plot a temporal slice IDT of the voxels image of tomobj X in pseudo
   Hounsfield unit ((µ-µ_water)/µ_water)*1000 on window NWIN.
 */
{
    x = 1000*((1.0/TRTT_WATER_ABSORPTION)*X.voxels - 1.0);
    mywindow, nwin; fma; img_plot, x(..,idt), cbar=1, vert=1, cmin=cmin, cmax=cmax, format="%.3f"; palette, "gray.gp";
}

func trtt_2D_plot_voxels_slice(x, idt, nwin, cmin=, cmax=, roi=)
/* DOCUMENT trtt_2D_plot_voxels_slice, x, idt, nwin, cmin=, cmax=

   Plot a temporal slice IDT of the voxels image X on window NWIN.
 */
{
    // xHU = 1000*((1.0/TRTT_WATER_ABSORPTION)*x - 1.0);
    if (!is_void(roi)) {
        mywindow, nwin; fma; img_plot, x(roi(1):roi(2),roi(3):roi(4),idt), cbar=1, vert=1, cmin=cmin, cmax=cmax;
    } else {
        mywindow, nwin; fma; img_plot, x(..,idt), cbar=1, vert=1, cmin=cmin, cmax=cmax;   
    }
    cmap, "gray";
}
trtt_plot_slice = trtt_2D_plot_voxels_slice;

func trtt_2D_plot_voxels(x, nwin, cmin=, cmax=, roi=)
/* DOCUMENT trtt_2D_plot_voxels, x, nwin, cmin=, cmax=

   Plot a static voxels image X on window NWIN.
 */
{
    // xHU = 1000*((1.0/TRTT_WATER_ABSORPTION)*x - 1.0);
    if (!is_void(roi)) {
        mywindow, nwin; fma; img_plot, x(roi(1):roi(2),roi(3):roi(4)), cbar=1, vert=1, cmin=cmin, cmax=cmax;
    } else {
        mywindow, nwin; fma; img_plot, x, cbar=1, vert=1, cmin=cmin, cmax=cmax;   
    }
    cmap, "gray";
}
trtt_plot_vox = trtt_2D_plot_voxels;

func trtt_2D_plot_voxels_2(x, nwin, cmin=, cmax=, roi=)
/* DOCUMENT trtt_2D_plot_voxels, x, nwin, cmin=, cmax=

   Plot a static voxels image X on window NWIN.
 */
{
    xHU = 1000*((1.0/TRTT_WATER_ABSORPTION)*x - 1.0);
    if (!is_void(roi)) {
        window, nwin, style="nobox.gs"; fma; img_plot, xHU(roi(1):roi(2),roi(3):roi(4)), cmin=cmin, cmax=cmax;
    } else {
        window, nwin, style="nobox.gs"; fma; img_plot, xHU, cmin=cmin, cmax=cmax;   
    }
    cmap, "gray";
}
trtt_plot_vox_2 = trtt_2D_plot_voxels_2;

func trtt_2D_plot_sinogram(Cops, nwin, cmin=, cmax=, dyn=)
/* DOCUMENT trtt_2D_plot_sinogram, Cops, nwin

   Plot sinogram data stored in the projectors set COPS on window NWIN.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    if (dyn)
        nv = h_get(Cops,Ck_list(1)).Ak.Yk.nv;
    else
        nv = h_get(Cops,Ck_list(1)).Yk.nv;
    sinogram=array(double,nv,ndata);

    for (k=1; k<=ndata; ++k) {
        if (dyn)
            Yk = h_get(Cops, Ck_list(k)).Ak.Yk;
        else
            Yk = h_get(Cops, Ck_list(k)).Yk;
        sinogram(,k) = Yk.pixels;
    }
    
    mywindow, nwin; fma; img_plot, sinogram, cbar=1, vert=1, cmin=cmin, cmax=cmax;
    cmap, "Blues";
}
trtt_plot_sino = trtt_2D_plot_sinogram;

func trtt_2D_plot_ref_sinogram(ya, nwin)
/* DOCUMENT trtt_2D_plot_ref_sinogram, ya, nwin

   Plot the NDATA reference projections YA on window NWIN.
 */
{
    for (k=1; k<=ndata; ++k) {
        mywindow, nwin; fma; img_plot, ya, cbar=1, vert=1, format="%.3f";
    }
}
trtt_plot_ref_sino = trtt_2D_plot_ref_sinogram;
