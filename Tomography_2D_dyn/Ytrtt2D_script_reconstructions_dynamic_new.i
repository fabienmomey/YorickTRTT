// mp_include, "Ytrtt2D_dyn.i";
include, "Ytrtt2D_dyn.i", 1;

/* TOMOBJ */
nx = 16;
ny = 16;
nt = 16;
xoff = 0.0; //FIXME: in mm
yoff = 0.0; //FIXME: in mm
s_scl = 1.0; //FIXME: in mm
s_deg = 3;
size_footprint = 10000;
/* TOMDATA */
nv = 512;
voff = 0.0; //FIXME: in mm
v_scl = 1.0; //FIXME: in mm

/* PROJECTOR */
mode = TRTT_FAN_BEAM;
type = TRTT_SPLINE_DRIVEN;
ndata = 300;
Rsc = 1000.0; //FIXME: in mm
// Rcd = 512.0; //FIXME: in mm
Rsd = 1536.0; // Rcd+Rsc;
SNR = 1.e3;

/* FIXME: THETA IS NOW THE ANGLE <Ox,OS> */
theta_step = 2*pi/ndata;
theta_range = theta_step*(ndata-1);
theta=trtt_span(0.0,theta_range,ndata);

/* Dynamic parameters */
tacq = 127.; //FIXME: in s
tcycl = 7.; //FIXME: in s
t_deg = 1;
t = (indgen(ndata)-1.0)*(tacq/(ndata-1.0));
t_index = (t%tcycl)*(nt/tcycl);
t_scl = tcycl/nt;

/* Noise features */                  
photon_flux = 2.e7; //FIXME: SNR² ~ photon_flux
read_noise = 5.;
add_noise = 1n;

/* ELLIPSES */
Nnorm = (nx-1)*s_scl*0.5;
// Ellipses = shepp_sparrow_semi_dynamic_2D(Nnorm);
Ellipses = shepp_sparrow_semi_dynamic_2D(Nnorm, tcycl=tcycl);
voxels_ref = trtt_create_dyn_phantom(Ellipses, nx, ny, s_scl, xoff, yoff, tacq, ndata);
voxels_refHU = 1000*((1.0/TRTT_WATER_ABSORPTION)*voxels_ref - 1.0);
/* Calculate analytic projection */
data = array(double, nv, ndata);
for (k=1; k<=ndata; ++k) {
    v = (indgen(nv)-0.5*(nv+1))*v_scl;
    v += voff;
    if (mode == TRTT_PARALLEL_BEAM) {
        data(,k) = R_ellipses_dyn_PB(Ellipses, v, nv, v_scl, -0.5*pi+theta(k), t(k), xoff, yoff);
    } else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
        data(,k) = R_ellipses_dyn_FB(Ellipses, v, nv, v_scl, -0.5*pi+theta(k), t(k), Rsc, Rsd, xoff, yoff);
    }
}

Htomo = trtt_2D_create_whole_simu_system(data,
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
                                         theta=theta,
                                         photon_flux=photon_flux,
                                         read_noise=read_noise,
                                         add_noise=add_noise,
                                         use_sparse_coefs=0n,
                                         ref_obj=voxels_ref);


local xref; eq_nocopy, xref, Htomo.Xref.voxels;
local Cops; eq_nocopy, Cops, Htomo.Cops;

mu_s=1.0;
mu_t=1.0;
regulTV = h_new(weight=[mu_s, mu_s, mu_t],
                threshold = 1.e-12,
                options = RGL_TOTVAR_ISOTROPIC);
XR = trtt_2D_optim_simu_launcher(Cops, trtt_cost_quadratic_opky, use_sparse_coefs=0n, x=array(double,nx,ny,nt), regulTV=regulTV, xmin=, mem=3, viewer=1n, win_viewer=4, win_viewer2=60, maxiter=100, maxeval=100, verbose=1n);

// mp_exec, "if (!mp_rank) quit;";
