// mp_include, "Ytrtt3D.i";
include, "Ytrtt3D.i", 1;

/* TOMOBJ */
nx = 16;
ny = 16;
nz = 16;
x_off = 0.0;
y_off = 0.0;
z_off = 0.0;
nt = 1;
s_scl = 1.0;
s_deg = 3;
size_footprint = 100000;

/* TOMDATA */
nu = 32;
nv = 32;
u_scl = 0.1;
v_scl = 0.1;
u_off = 0.0;
v_off = 0.0;

/* PROJECTOR */
mode = TRTT_CONE_BEAM;
type = TRTT_SPLINE_DRIVEN;
ndata = 720;
Rsc = 512.;
Rcd = 512.;
Rsd = Rcd+Rsc;

/* Noise features */                  
// photon_flux = 2.e7; //FIXME: SNR² ~ photon_flux
// read_noise = 5.;
// sig_noise_uniform=0.1*0.025*TRTT_WATER_ABSORPTION;
// var_noise_uniform=sig_noise_uniform^2.0;
// add_noise = 1n;

/* FIXME: THETA IS NOW THE ANGLE <Ox,OS> */
theta_step = 2*pi/ndata;
theta_range = theta_step*(ndata-1);
theta=trtt_span(0.0,theta_range,ndata);

/* ELLIPSES */
// Ellipses = shepp_logan_3D();
data_file = "Simul_Data/proj_boule.data";
object_file = "Simul_Data/boule.data";
r_boule = 0.5;
each = 720/ndata;
data = double(raw_read(data_file, float, nu, nv, 720, offset=512))(,,::each)(..,1:ndata);
voxels_ref = double(raw_read(object_file, float, nx, ny, nz, offset=512));
for(k=1; k<=ndata; ++k) {
    data(..,k)=transpose(data(::-1,,k));
}


Htomo = trtt_3D_create_whole_system(data,
                                    nu,
                                    nv,
                                    u_scl,
                                    v_scl,
                                    nx,
                                    ny,
                                    nz,
                                    s_scl,
                                    s_deg,
                                    size_footprint,
                                    ndata,
                                    mode,
                                    type,
                                    Rsc,
                                    Rsd,
                                    x_off=x_off,
                                    y_off=y_off,
                                    z_off=z_off,
                                    u_off=u_off,
                                    v_off=v_off,
                                    theta=theta,
                                    ref_obj=voxels_ref);

local xref; eq_nocopy, xref, Htomo.Xref.voxels;
local Cops; eq_nocopy, Cops, Htomo.Cops;

trtt3D_plot_sino, Cops, 0, wait=10;
data=trtt3D_get_sino(Cops);
// trtt3D_plot_vox, xref, 3, 1, cmin=0, cmax=30, wait=100;
mywindow, 1; fma; img_plot, xref(..,8), cmin=0, cmax=0.005, cbar=1, vert=1;

data_alg = trtt3D_get_ref_proj(Htomo);

/* GO RECONSTRUCTION ! */
// mu =(1./3.)*sig_noise_uniform*r_boule; regulTV = h_new(weight=[mu, mu, mu], threshold = 1.e-12, options = RGL_TOTVAR_ISOTROPIC, lweight = array(1.0,nx,ny,nz)); XR = trtt_3D_optim_simu_launcher(Cops, trtt3D_cost_quadratic_mpy, regulTV=regulTV, maxiter=1000, verbose=1);
// trtt3D_plot_vox, XR.x, 3, 3, cmin=0.0, cmax=30.0, wait=50;
// mywindow, 3; fma; img_plot, XR.x(..,8), cmin=0.0, cmax=0.005, cbar=1, vert=1;

// trtt_2D_save, Htomo, "", overwrite=1;
