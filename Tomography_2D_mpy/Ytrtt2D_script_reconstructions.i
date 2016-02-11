mp_include, "Ytrtt2D_mpy.i", 1;

/* TOMOBJ */
nx = 256;
ny = 256;
xoff = 0; //-100; //FIXME: in pix
yoff = 0; //FIXME: in pix
s_scl = 0.1; //FIXME: in cm
s_deg = 3;
size_footprint = 10000;

/* TOMDATA */
nv = 512;
voff = 0.0; //FIXME: in cm
v_scl = 0.1; //FIXME: in cm
t_index = 0.0;

/* PROJECTOR */
mode = TRTT_FAN_BEAM;
type = TRTT_SPLINE_DRIVEN;
ndata = 60;
Rsc = 100.0; //FIXME: in cm
// Rcd = 51.2; //FIXME: in cm
Rsd = 153.6; // Rcd+Rsc;
SNR = 1.e3;

/* Noise features */                  
photon_flux = 1.e8; //FIXME: SNR² ~ photon_flux
read_noise = 5.;

/* FIXME: THETA IS NOW THE ANGLE <Ox,OS> */
theta_step = 2*pi/ndata;
theta_range = theta_step*(ndata-1);
theta=trtt_span(0.0,theta_range,ndata);

/* ELLIPSES */
Nnorm = (nx-1)*s_scl*0.5;
Ellipses = shepp_sparrow_2D(Nnorm);
voxels_ref = trtt_create_ellipses_phantom(Ellipses, nx, ny, s_scl, xoff, yoff);
/* Calculate analytic projection */
data = array(double, nv, ndata);
for (k=1; k<=ndata; ++k) {
    v = (indgen(nv)-0.5*(nv+1))*v_scl;
    v += voff;
    if (mode == TRTT_PARALLEL_BEAM) {
        data(,k) = R_ellipses_PB(Ellipses, v, nv, v_scl, -0.5*pi+theta(k), s_scl, xoff, yoff);
    } else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
        data(,k) = R_ellipses_FB(Ellipses, v, nv, v_scl, -0.5*pi+theta(k), Rsc, Rsd, s_scl, xoff, yoff);
    }
}

Htomo = trtt_2D_create_whole_simu_system(data,
                                         nv,
                                         v_scl,
                                         400,
                                         400,
                                         s_scl,
                                         s_deg,
                                         size_footprint,
                                         ndata,
                                         mode,
                                         type,
                                         Rsc,
                                         Rsd,
                                         // xoff=xoff,
                                         // yoff=yoff,
                                         // voff=voff,
                                         theta=theta,
                                         photon_flux=photon_flux,
                                         read_noise=read_noise,
                                         add_noise=add_noise,
                                         ref_obj=voxels_ref,
                                         matrix=0);

xref = Htomo.Xref.voxels;
Cops = Htomo.Cops;
// sino_analytic = trtt_get_sino(Cops);
// sino_algebraic = trtt_get_ref_sino(Htomo);
/* SHEPP LOGAN SIMULE */
// trtt_plot_vox, Htomo.Xref.voxels, 1, 0, cmin=-100, cmax=100;
// trtt_plot_sino, Htomo.Cops, 1;
// trtt_plot_ref_sino, sino_algebraic, 2;


/* GO RECONSTRUCTION ! */
mu = 500;
regulTV = h_new(weight=[mu, mu], threshold = 1.e-5, options = RGL_TOTVAR_ISOTROPIC);
XR = trtt_2D_optim_simu_launcher(Cops, trtt_cost_quadratic_mpy_opky, use_sparse_coefs=0n, regulTV=regulTV, viewer=1n, win_viewer=2, win_viewer2=60, cmin=, cmax=, maxiter=500, verbose=1); trtt_plot_vox, XR.x, 3, cmin=-100, cmax=100; palette, "gray.gp";
X


