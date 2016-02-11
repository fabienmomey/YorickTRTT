mp_include, "Ytrtt2D_mpy.i";

write,format="MP_SIZE = %d \n", mp_size;

data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320853170.4/";
DATA = yhd_restore(data_dir+"Data4TRTT");

ntheta = DATA.ndata;
Rsc = DATA.Rsc;
Rsd = DATA.Rsd;
theta = DATA.theta/*(5:0)*/; //FIXME: the first five projections are static
ntheta = numberof(theta);
local data; eq_nocopy, data, DATA.data/*(..,5:)*/;

// for (k=1; k<=ntheta; ++k) {
//     mywindow, 0; fma; img_plot, data(..,k), cbar=1, vert=1, format="%.3f";
//     pause, 1;
// }

/* TOMOBJ */
nx = 190;
ny = 190;
xoff = -45.0; //FIXME: in mm
yoff = -20.0;   //FIXME: in mm
s_scl = 2.0; //FIXME: in mm
s_deg = 3;
size_footprint = 10000;

/* TOMDATA */
each = 1; //FIXME: we take ~ 1 projection of EACH
nv = DATA.nv;
voff = DATA.voff/*(5:0)*/; //FIXME: in mm
uoff = DATA.uoff/*(5:0)*/;
v_scl = DATA.v_scl; //FIXME: in mm
t_index = 0.0;

uoff = long(floor((uoff/v_scl) +0.5));

/* PROJECTOR */
mode = TRTT_FAN_BEAM;
type = TRTT_SPLINE_DRIVEN;

/* Noise features */                  
photon_flux = 1.e6; //FIXME: SNR² ~ photon_flux
read_noise = 5.;
var_noise_uniform = 1.;

/* Extract angles */
theta = theta(::each);
voff = voff(::each);
uoff = uoff(::each);
ndata = numberof(theta);

/* Extract data */
data2d = array(double,nv,ndata);
id_theta = indgen(ntheta)(::each);
midpanel = long(floor(nv/2));
for (k=1; k<=ndata; ++k) {
    data2d(,k) = data(midpanel-uoff(k),,id_theta(k));
}

id1 = where(theta>pi);
id2 = where(theta<=pi);
theta(id2) = -1.0*theta(id2);
theta(id1) = 2*pi-1.0*theta(id1);

DATA=[];
data=[];

Htomo = trtt_2D_create_whole_simu_system(data2d,
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
                                         xoff=xoff,
                                         yoff=yoff,
                                         voff=-voff, /*FIXME:FIXME: IL FAUT METTRE -VOFF !!!*/
                                         theta=theta,
                                         // var_noise_uniform=var_noise_uniform,
                                         // photon_flux=photon_flux,
                                         // read_noise=read_noise,
                                         use_sparse_coefs=0n);


local Cops; eq_nocopy, Cops, Htomo.Cops;
Ck_list=Cops.Ck_list;

trtt_plot_sino, Htomo.Cops, 1;

/* GO RECONSTRUCTION ! */
dweights = array(double,nv,ndata);
dweights(1:511,)=1.0;
id_invalid = where(data2d == -1.0);
if (is_array(id_invalid))
    dweights(id_invalid) = 0.0;

mu =1.0;
regulTV = h_new(weight=mu,
                threshold = 1.e-2,
                options = RGL_TOTVAR_ISOTROPIC);

XR = trtt_2D_optim_launcher(Cops, trtt_cost_quadratic_mpy_opky, x=array(double,nx,ny), use_sparse_coefs=0n, dweights=dweights, mem=5, regulTV=regulTV, xmin=0.0, viewer=1n, win_viewer=4, win_viewer2=60, maxiter=20, maxeval=20, verbose=1n);
trtt_plot_vox, XR.x, 5;
