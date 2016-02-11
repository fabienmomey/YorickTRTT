/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt3D_plot.i --
 *
 * TRTT 3D plot function package.
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

/* USER FUNCTIONS ============================================================ */

func trtt_3D_get_sinogram(Cops, dyn=)
/* DOCUMENT trtt_3D_get_sinogram, Cops, dyn=

   Get sinogram data stored in the projectors set COPS, and arrange it in a
   single multidimensional array of dimensions (nu x nv x ndata).
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    if (!dyn) {
        nu = h_get(Cops,Ck_list(1)).Yk.nu;
        nv = h_get(Cops,Ck_list(1)).Yk.nv;
    } else {
        nu = h_get(Cops,Ck_list(1)).Ak.Yk.nu;
        nv = h_get(Cops,Ck_list(1)).Ak.Yk.nv;
    }
    sinogram=array(double,nu,nv,ndata);

    for (k=1; k<=ndata; ++k) {
        if (!dyn) {
            Yk = h_get(Cops, Ck_list(k)).Yk;
        } else {
            Yk = h_get(Cops, Ck_list(k)).Ak.Yk;
        }
        sinogram(,,k) = Yk.pixels;
    }

    return sinogram;
}
trtt3D_get_sino = trtt_3D_get_sinogram;

func trtt_3D_get_ref_projections(Htomo)
/* DOCUMENT trtt_3D_get_ref_projections, Htomo

   Get projections of the reference voxel image XREF "passed through" the
   projectors Ck contained in COPS.
 */
{
    X = Htomo.X;
    Y = Htomo.Y;
    Xref = Htomo.Xref;
    Cops = Htomo.Cops;
    s_deg = X.s_deg;
    t_deg = X.t_deg;
    x = Xref.voxels;
    dims=dimsof(x);
    if (type == TRTT_SPLINE_DRIVEN) {
        if (dims(1)==4)
            c = spl_spline_samples_to_coefficients(x, s_deg, s_deg, s_deg, t_deg);
        else
            c = spl_spline_samples_to_coefficients(x, s_deg, s_deg, s_deg);
    }
    Ck_list = Cops.Ck_list;
    Yk_list = Y.Yk_list;
    ndata = numberof(Ck_list);
    nu = h_get(Y,Yk_list(1)).nu;
    nv = h_get(Y,Yk_list(1)).nv;
    ya = array(double,nu,nv,ndata);
    for (k=1; k<=ndata; ++k) {
        write, format="Projection n°%d\n", k;
        Ck = h_get(Cops,Ck_list(k));
        Yk = h_get(Y,Yk_list(k));
        if (type == TRTT_SPLINE_DRIVEN)
            ya(..,k) = Ck(c);
        else
            ya(..,k) = Ck(x);
    }

    return ya;
}
trtt3D_get_ref_proj = trtt_3D_get_ref_projections;

func trtt_3D_plot_sinogram(Cops, nwin, wait=)
/* DOCUMENT trtt_3D_plot_sinogram, Cops, nwin, wait=

   Plot projections stored in the projectors set COPS on window NWIN. WAIT
   specify the pause (in ms) between each display.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);

    if (is_void(wait)) wait=500;

    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Cops, Ck_list(k)).Yk;
        mywindow, nwin; fma; img_plot, Yk.pixels, cbar=1, vert=1, format="%.3f";
        pause, wait;
    }
}
trtt3D_plot_sino = trtt_3D_plot_sinogram;

func trtt_3D_plot_ref_projections(proj, ndata, nwin, wait=)
/* DOCUMENT trtt_3D_plot_ref_projections, proj, ndata, nwin, wait=

   Plot the NDATA reference projections PROJ on window NWIN. WAIT specify
   the pause (in ms) between each display.
 */
{
    if (is_void(wait)) wait=500;

    for (k=1; k<=ndata; ++k) {
        mywindow, nwin; fma; img_plot, proj(..,k), cbar=1, vert=1, format="%.3f";
        pause, wait;
    }
}
trtt3D_plot_ref_proj = trtt_3D_plot_ref_projections;

func trtt_3D_plot_voxels(x, axis, k, nwin, cmin=, cmax=)
/* DOCUMENT trtt_3D_plot_voxels, x, axis, nwin, cmin=, cmax=, wait=

   Plot slices of object voxels X, according to the varying axis AXIS:
   - 1 for x axis.
   - 2 for y axis.
   - 3 for z axis.
   The plot is onde on window NWIN. CMIN and CMAX are value bounds. WAIT specify
   the pause (in ms) between each display.

   SEE ALSO: img_plot.
 */
{
    // xHU = 1000*((1.0/TRTT_WATER_ABSORPTION)*x - 1.0);
    // dims=dimsof(x);
    // nx=dims(2); ny=dims(3); nz=dims(4);

    // if (is_void(wait)) wait=100;

    /* Varying X axis */
    if (axis == 1) {
        // for (k=1; k<=nx; ++k) {
            mywindow, nwin; fma; img_plot, x(k,,), cbar=1, vert=1, format="%.3f", cmin=cmin, cmax=cmax;
            // pause, wait;
        // }
    }
    /* Varying Y axis */
    if (axis == 2) {
        // for (k=1; k<=ny; ++k) {
            mywindow, nwin; fma; img_plot, x(,k,), cbar=1, vert=1, format="%.3f", cmin=cmin, cmax=cmax;
            // pause, wait;
        // }
    }
    /* Varying Z axis */
    if (axis == 3) {
        // for (k=1; k<=nz; ++k) {
            mywindow, nwin; fma; img_plot, x(,,k), cbar=1, vert=1, format="%.3f", cmin=cmin, cmax=cmax;
            // pause, wait;
        // }
    }
    cmap, "gray";
}
trtt3D_plot_vox = trtt_3D_plot_voxels;

func trtt_3D_plot_slice(x, axis, k, idt, nwin, cmin=, cmax=)
/* DOCUMENT trtt_3D_plot_voxels, voxels, k, axis, nwin, cmin=, cmax=

   Plot a particular slice K of the VOXELS, according to the varying axis AXIS:
   - 1 for x axis.
   - 2 for y axis.
   - 3 for z axis.
   The plot is onde on window NWIN. CMIN and CMAX are value bounds.

   SEE ALSO: img_plot.
 */
{
    // xHU = 1000*((1.0/TRTT_WATER_ABSORPTION)*voxels - 1.0);
    /* Varying X axis */
    if (axis == 1) {
        mywindow, nwin; fma; img_plot, x(k,,,idt), cbar=1, vert=1, format="%.3f", cmin=cmin, cmax=cmax;
    }
    /* Varying Y axis */
    if (axis == 2) {      
        mywindow, nwin; fma; img_plot, x(,k,,idt), cbar=1, vert=1, format="%.3f", cmin=cmin, cmax=cmax;
    }
    /* Varying Z axis */
    if (axis == 3) {     
        mywindow, nwin; fma; img_plot, x(,,k,idt), cbar=1, vert=1, format="%.3f", cmin=cmin, cmax=cmax;
    }
    cmap, "gray";
}
trtt3D_plot_slice = trtt_3D_plot_slice;
