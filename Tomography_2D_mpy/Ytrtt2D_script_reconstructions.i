// mp_include, "Ytrtt2D_mpy.i";
include, "../Tomography_2D/Ytrtt2D.i", 1;

/* TOMOBJ */
nx = 300;
ny = 300;
xoff = 0; //-100; //FIXME: in pix
yoff = 0; //FIXME: in pix
s_scl = 256./300; //FIXME: in mm
s_deg = 0;
size_footprint = 10000;

/* TOMDATA */
nv = 512;
voff = 0.0; //FIXME: in mm
v_scl = 1.0; //FIXME: in mm
t_index = 0.0;

/* PROJECTOR */
mode = TRTT_FAN_BEAM;
type = TRTT_SPLINE_DRIVEN;
ndata = 300;
Rsc = 1000.; //FIXME: in mm
// Rcd = 51.2; //FIXME: in mm
Rsd = 1536.; // Rcd+Rsc;
SNR = 1.e3;

/* Noise features */                  
photon_flux = 1.e8; //FIXME: SNR² ~ photon_flux
read_noise = 5.;

/* FIXME: THETA IS NOW THE ANGLE <Ox,OS> */
theta_step = 2*pi/ndata;
theta_range = theta_step*(ndata-1);
theta=trtt_span(0.0,theta_range,ndata);
// theta+⁼theta_step/3.0;

/* ELLIPSES */
Nnorm = (nx-1)*s_scl*0.5;
Ellipses = shepp_logan_2D(Nnorm);
voxels_ref = trtt_create_ellipses_phantom(Ellipses, nx, ny, s_scl, xoff, yoff);
/* Calculate analytic projection */
data = array(double, nv, ndata);
for (k=1; k<=ndata; ++k) {
    v = (indgen(nv)-0.5*(nv+1))*v_scl;
    v += voff;
    if (mode == TRTT_PARALLEL_BEAM) {
        data(,k) = R_ellipses_PB(Ellipses, v, nv, v_scl, -0.5*pi+theta(k), xoff, yoff);
    } else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
        data(,k) = R_ellipses_FB(Ellipses, v, nv, v_scl, -0.5*pi+theta(k), Rsc, Rsd, xoff, yoff);
    }
}

Htomo = trtt_2D_create_whole_simu_system(data,
                                         nv,
                                         v_scl,
                                         nx,
                                         ny,
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
                                         // photon_flux=photon_flux,
                                         // read_noise=read_noise,
                                         // add_noise=add_noise,
                                         ref_obj=voxels_ref,
                                         use_sparse_coefs=0n);

xref = Htomo.Xref.voxels;
Cops = Htomo.Cops;
sino_analytic = trtt_get_sino(Cops);
sino_algebraic = trtt_get_ref_sino(Htomo);
/* SHEPP LOGAN SIMULE */
// trtt_plot_vox, Htomo.Xref.voxels, 0;
// trtt_plot_sino, Htomo.Cops, 1;
// trtt_plot_ref_sino, sino_algebraic, 2;


/* GO RECONSTRUCTION ! */
// mu = 0.1;
// regulTV = h_new(weight=[mu, mu], threshold = 1.e-3, options = RGL_TOTVAR_ISOTROPIC);
// XR = trtt_2D_optim_simu_launcher(Cops, trtt_cost_quadratic_mpy_opky, use_sparse_coefs=0n, mem=5, regulTV=, viewer=1n, win_viewer=2, win_viewer2=60, cmin=1.0, cmax=1.05, maxiter=20, verbose=1); trtt_plot_vox, XR.x, 3, cmin=1.0, cmax=1.05; palette, "gray.gp";



